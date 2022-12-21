# This script performs the following tasks:
# 1) download weather stations with more than 100 years of precipitation or temperature data
# 2) calculate the residuals from each station using ctbi (no outliers flagged) with 3 bin intervals : 5,10 and 20 days
# 3) calculate the % of points removed for each method on all available points per station
# 4) create the histograms by pooling all anomalies from all stations (for vizualization purpose only)
# 5) calculate the % of points removed for each method for sample sizes varying from 10 to 100


path.main <- "~/Science/ctbi_temp2/"
setwd(path.main)

set.seed(1)

name.prcp.out <- 'thresholds_precipitation_v5.csv'
name.temp.out <- 'thresholds_temperature_v5.csv'
name.hist.prcp.out <- 'hist_precipitation_v5.csv'
name.hist.temp.out <- 'hist_temperature_v5.csv'
name.small.prcp.out <- 'thresholds_precipitation_smallsample_v7.csv'
name.small.temp.out <- 'thresholds_temperature_smallsample_v7.csv'


setwd('C:/Users/fritter/Documents/GitHub/ctbi_article/PartI.outliers/step4_comparisonofmethods/')
source('f.Kimber_v2.R')
source('f.Hubert_v2.R')
source('f.LogBox_v17.R')
source('f.Leys_v2.R')
source('f.Schwertman_v2.R')
source('f.Barbato_v1.R')
setwd(path.main)

Z.schwertman <- 3
if(1)
{
  # need to find the k.n related to n=128 and n=256 for Schwertman with a linear interpolation
  x1 <- 100
  x2 <- 200
  y1 <- 1.34588
  y2 <- 1.34740
  slope0 <- (y2-y1)/(x2-x1)
  intercept0 <- y1 - slope0*x1
  y.128 <- slope0*128+intercept0
  x1 <- 200
  x2 <- 300
  y1 <- 1.34740
  y2 <- 1.34792
  slope0 <- (y2-y1)/(x2-x1)
  intercept0 <- y1 - slope0*x1
  y.256 <- slope0*256+intercept0

  # constants used in the Schwertman model to adjust for sample size effect
  k.n.Schwertman <- c(1.33318,1.34004,1.34424,y.128,y.256,1.34898,1.34898,1.34898,1.34898,1.34898,1.34898)
}

library(data.table)
library(rnoaa)
library(ctbi)
library(mrfDepth) # for the med couple

# get all stations
if(0)
{
  if(1) # daily precipitation
  {
    data <- fread('GHCN_stationlist.csv')
    data <- data[element %in% c('PRCP'),]
    data[,diff_year := (last_year-first_year)]
    data <- data[diff_year >= 100,] # select stations with at least 100 years of data
    list.stations <- as.character(unlist(unique(data[,id])))
    n.station <- length(list.stations)
    for(i in 1:n.station)
    {
      print(i)
      df.tp <- ghcnd_search(stationid = list.stations[i],var = 'PRCP',refresh = FALSE)
      w.loop <- setDT(df.tp$prcp)

      if(sum(names(w.loop) %in% c('prcp','qflag')) == 2)
      {
        w.loop[,prcp := round(prcp/10,digits=2)]
        w.loop[qflag != ' ',prcp := NA]
        w.loop <- w.loop[,c('date','prcp'),with=F]

        if(sum(!is.na(w.loop[,prcp])) > 10000)
        {
          fwrite(w.loop,file= paste0(path.main,'IN_precipitation/',list.stations[i],'.csv'))
        }
      }
    }
  }

  if(1) # daily temperature
  {
    data <- fread('GHCN_stationlist.csv')
    data <- data[element %in% c('TMIN','TMAX'),]
    data[,diff_year := (last_year-first_year)]
    data <- data[diff_year >= 100,] # select stations with at least 100 years of data
    list.stations <- as.character(unlist(unique(data[,id])))
    n.station <- length(list.stations)
    for(i in 1:n.station)
    {
      print(i)
      df.tp <- ghcnd_search(stationid = list.stations[i],var = c('TMIN','TMAX'),refresh = FALSE)

      df.min <- setDT(df.tp$tmin)
      df.max <- setDT(df.tp$tmax)

      logmin <- sum(names(df.min) %in% c('tmin','qflag')) == 2
      logmax <- sum(names(df.max) %in% c('tmax','qflag')) == 2

      if(logmin & logmax)
      {
        df.min[qflag != ' ',tmin := NA]
        df.max[qflag != ' ',tmax := NA]

        w.loop <- merge(df.min[,c('date','tmin'),with=FALSE],df.max[,c('date','tmax'),with=FALSE],by='date')
        w.loop[,tavg := round((tmin+tmax)/20,digits=2)]
        w.loop[,':='(tmin=NULL,tmax=NULL)]

        if(sum(!is.na(w.loop[,tavg])) > 10000)
        {
          fwrite(w.loop,file= paste0(path.main,'IN_temperature/',list.stations[i],'.csv'))
        }
      }
    }
  }
}

# calculate the precipitation & temperature residuals for three different long.term fits (bins with 5,10 and 20 points)
if(0)
{
  if(1) # precipitation residuals
  {
    setwd(paste0(path.main,'IN_precipitation/'))

    list.stations <- list.files(pattern = "\\.csv$")
    n.station <- length(list.stations)

    for(i in 1:n.station)
    {
      print(c(i,n.station))

      y.raw <- fread(list.stations[i])

      list.main5 <- ctbi(copy(y.raw),bin.period='5 days',bin.side=as.Date('2000-01-01'),coeff.outlier = NA,SCI.min = NA)
      y <- list.main5$data0[,c('date','prcp','index.bin','long.term','cycle'),with=F]
      y <- y[index.bin > 0,]
      y <- y[!is.na(prcp),]
      y <- y[prcp > 0,]
      y <- y[,resid := (prcp-long.term-cycle)]
      y5 <- unlist(y[!is.na(resid),resid],use.names = F)

      list.main5 <- ctbi(copy(y.raw),bin.period='10 days',bin.side=as.Date('2000-01-01'),coeff.outlier = NA,SCI.min = NA)
      y <- list.main5$data0[,c('date','prcp','index.bin','long.term','cycle'),with=F]
      y <- y[index.bin > 0,]
      y <- y[!is.na(prcp),]
      y <- y[prcp > 0,]
      y <- y[,resid := (prcp-long.term-cycle)]
      y10 <- unlist(y[!is.na(resid),resid],use.names = F)

      list.main5 <- ctbi(copy(y.raw),bin.period='20 days',bin.side=as.Date('2000-01-01'),coeff.outlier = NA,SCI.min = NA)
      y <- list.main5$data0[,c('date','prcp','index.bin','long.term','cycle'),with=F]
      y <- y[index.bin > 0,]
      y <- y[!is.na(prcp),]
      y <- y[prcp > 0,]
      y <- y[,resid := (prcp-long.term-cycle)]
      y20 <- unlist(y[!is.na(resid),resid],use.names = F)

      n5 <- length(y5)
      n10 <- length(y10)
      n20 <- length(y20)

      n.max <- max(n5,n10,n20)

      if(n.max > 16)
      {
        data.out <- data.frame(y5=rep(NA,n.max),y10=rep(NA,n.max),y20=rep(NA,n.max),stringsAsFactors = F)
        data.out[1:n5,'y5'] <- round(y5,digits=4)
        data.out[1:n10,'y10'] <- round(y10,digits=4)
        data.out[1:n20,'y20'] <- round(y20,digits=4)

        fwrite(data.out,file=paste0(path.main,'/IN_precipitation_residuals/',list.stations[i]),sep=',')
      }
    }
  }

  if(1) # temperature residuals
  {
    setwd(paste0(path.main,'IN_temperature/'))

    list.stations <- list.files(pattern = "\\.csv$")
    n.station <- length(list.stations)

    for(i in 1:n.station)
    {
      print(c(i,n.station))

      y.raw <- fread(list.stations[i])

      list.main5 <- ctbi(copy(y.raw),bin.period='5 days',bin.side=as.Date('2000-01-01'),coeff.outlier = NA,SCI.min = NA)
      y <- list.main5$data0[,c('date','tavg','index.bin','long.term','cycle'),with=F]
      y <- y[index.bin > 0,]
      y <- y[!is.na(tavg),]
      y <- y[,resid := (tavg-long.term-cycle)]
      y5 <- unlist(y[!is.na(resid),resid],use.names = F)

      list.main5 <- ctbi(copy(y.raw),bin.period='10 days',bin.side=as.Date('2000-01-01'),coeff.outlier = NA,SCI.min = NA)
      y <- list.main5$data0[,c('date','tavg','index.bin','long.term','cycle'),with=F]
      y <- y[index.bin > 0,]
      y <- y[!is.na(tavg),]
      y <- y[,resid := (tavg-long.term-cycle)]
      y10 <- unlist(y[!is.na(resid),resid],use.names = F)

      list.main5 <- ctbi(copy(y.raw),bin.period='20 days',bin.side=as.Date('2000-01-01'),coeff.outlier = NA,SCI.min = NA)
      y <- list.main5$data0[,c('date','tavg','index.bin','long.term','cycle'),with=F]
      y <- y[index.bin > 0,]
      y <- y[!is.na(tavg),]
      y <- y[,resid := (tavg-long.term-cycle)]
      y20 <- unlist(y[!is.na(resid),resid],use.names = F)

      n5 <- length(y5)
      n10 <- length(y10)
      n20 <- length(y20)

      n.max <- max(n5,n10,n20)

      if(n.max > 16)
      {
        data.out <- data.frame(y5=rep(NA,n.max),y10=rep(NA,n.max),y20=rep(NA,n.max),stringsAsFactors = F)
        data.out[1:n5,'y5'] <- round(y5,digits=4)
        data.out[1:n10,'y10'] <- round(y10,digits=4)
        data.out[1:n20,'y20'] <- round(y20,digits=4)

        fwrite(data.out,file=paste0(path.main,'/IN_temperature_residuals/',list.stations[i]),sep=',')
      }
    }
  }
}

# calculate the number of points removed for each station (all available points)
if(1)
{
  if(1) # precipitation residuals
  {
    setwd(paste0(path.main,'IN_precipitation_residuals/'))
    list.stations <- list.files(pattern = "\\.csv$")
    n.station <- length(list.stations)
    dtf.m.P <- data.frame(stations=substr(list.stations,1,11),n=rep(NA,n.station),kim=rep(NA,n.station),hub=rep(NA,n.station),sch=rep(NA,n.station),ley=rep(NA,n.station),logbox=rep(NA,n.station),bar=rep(NA,n.station),l.kim=rep(NA,n.station),l.hub=rep(NA,n.station),l.sch=rep(NA,n.station),l.ley=rep(NA,n.station),l.logbox=rep(NA,n.station),l.bar=rep(NA,n.station),u.kim=rep(NA,n.station),u.hub=rep(NA,n.station),u.sch=rep(NA,n.station),u.ley=rep(NA,n.station),u.logbox=rep(NA,n.station),u.bar=rep(NA,n.station),stringsAsFactors = F)
    for(i in 1:n.station)
    {
      print(c(i,n.station))
      d <- read.table(list.stations[i],sep=',',header=T,stringsAsFactors = F)
      y5 <- d[,1]
      y10 <- d[,2]
      y20 <- d[,3]
      y5 <- y5[!is.na(y5)]
      y10 <- y10[!is.na(y10)]
      y20 <- y20[!is.na(y20)]

      n5 <- length(y5)
      n10 <- length(y10)
      n20 <- length(y20)

      n.all <- round(mean(c(n5,n10,n20)),digits=4)
      n.out.kim <- c()
      n.out.hub <- c()
      n.out.sch <- c()
      n.out.ley <- c()
      n.out.logbox <- c()
      n.out.bar <- c()
      
      l.kim <- c()
      l.hub <- c()
      l.sch <- c()
      l.ley <- c()
      l.logbox <- c()
      l.bar <- c()
      
      u.kim <- c()
      u.hub <- c()
      u.sch <- c()
      u.ley <- c()
      u.logbox <- c()
      u.bar <- c()
      
      for(j in 1:3)
      {
        if(j==1)
        {
          d.loop <- y5
          n.loop <- n5
        }
        if(j==2)
        {
          d.loop <- y10
          n.loop <- n10
        }
        if(j==3)
        {
          d.loop <- y20
          n.loop <- n20
        }

        ###
        MC <- as.numeric(medcouple(d.loop))
        MAD0 <- as.numeric(mad(d.loop))
        q0.125 <- as.numeric(quantile(d.loop,0.125))
        q0.25 <- as.numeric(quantile(d.loop,0.25))
        q0.375 <- as.numeric(quantile(d.loop,0.375))
        q0.50 <- as.numeric(quantile(d.loop,0.50))
        q0.625 <- as.numeric(quantile(d.loop,0.625))
        q0.75 <- as.numeric(quantile(d.loop,0.75))
        q0.875 <- as.numeric(quantile(d.loop,0.875))
        coeff.outlier <- 'auto'
        
        nlu.kim <- f.Kimber(d.loop,q0.25,q0.50,q0.75)
        nlu.hub <- f.Hubert(d.loop,q0.25,q0.75,MC)
        nlu.sch <- f.Schwertman(d.loop,q0.25,q0.50,q0.75,Z.schwertman,k.n.Schwertman[11])
        nlu.ley <- f.Leys(d.loop,q0.50,MAD0)
        nlu.logbox <- f.LogBox(d.loop,coeff.outlier,q0.125,q0.25,q0.375,q0.625,q0.75,q0.875)
        nlu.bar <- f.Barbato(d.loop,q0.25,q0.75)
        

        n.out.kim <- c(n.out.kim,nlu.kim[1])
        n.out.hub <- c(n.out.hub,nlu.hub[1])
        n.out.sch <- c(n.out.sch,nlu.sch[1])
        n.out.ley <- c(n.out.ley,nlu.ley[1])
        n.out.logbox <- c(n.out.logbox,nlu.logbox[1])
        n.out.bar <- c(n.out.bar,nlu.bar[1])
        
        l.kim <- c(l.kim,nlu.kim[2])
        l.hub <- c(l.hub,nlu.hub[2])
        l.sch <- c(l.sch,nlu.sch[2])
        l.ley <- c(l.ley,nlu.ley[2])
        l.logbox <- c(l.logbox,nlu.logbox[2])
        l.bar <- c(l.bar,nlu.bar[2])
        
        u.kim <- c(u.kim,nlu.kim[3])
        u.hub <- c(u.hub,nlu.hub[3])
        u.sch <- c(u.sch,nlu.sch[3])
        u.ley <- c(u.ley,nlu.ley[3])
        u.logbox <- c(u.logbox,nlu.logbox[3])
        u.bar <- c(u.bar,nlu.bar[3])
      }

      n.out.kim <- round(mean(n.out.kim),digits=4)
      n.out.hub <- round(mean(n.out.hub),digits=4)
      n.out.sch <- round(mean(n.out.sch),digits=4)
      n.out.ley <- round(mean(n.out.ley),digits=4)
      n.out.logbox <- round(mean(n.out.logbox),digits=4)
      n.out.bar <- round(mean(n.out.bar),digits=4)

      l.kim <- round(mean(l.kim),digits=4)
      l.hub <- round(mean(l.hub),digits=4)
      l.sch <- round(mean(l.sch),digits=4)
      l.ley <- round(mean(l.ley),digits=4)
      l.logbox <- round(mean(l.logbox),digits=4)
      l.bar <- round(mean(l.bar),digits=4)
      
      u.kim <- round(mean(u.kim),digits=4)
      u.hub <- round(mean(u.hub),digits=4)
      u.sch <- round(mean(u.sch),digits=4)
      u.ley <- round(mean(u.ley),digits=4)
      u.logbox <- round(mean(u.logbox),digits=4)
      u.bar <- round(mean(u.bar),digits=4)
      
      dtf.m.P[i,'n'] <- n.all
      dtf.m.P[i,'kim'] <- n.out.kim
      dtf.m.P[i,'hub'] <- n.out.hub
      dtf.m.P[i,'sch'] <- n.out.sch
      dtf.m.P[i,'ley'] <- n.out.ley
      dtf.m.P[i,'logbox'] <- n.out.logbox
      dtf.m.P[i,'bar'] <- n.out.bar
      
      dtf.m.P[i,'l.kim'] <- l.kim
      dtf.m.P[i,'l.hub'] <- l.hub
      dtf.m.P[i,'l.sch'] <- l.sch
      dtf.m.P[i,'l.ley'] <- l.ley
      dtf.m.P[i,'l.logbox'] <- l.logbox
      dtf.m.P[i,'l.bar'] <- l.bar
      
      dtf.m.P[i,'u.kim'] <- u.kim
      dtf.m.P[i,'u.hub'] <- u.hub
      dtf.m.P[i,'u.sch'] <- u.sch
      dtf.m.P[i,'u.ley'] <- u.ley
      dtf.m.P[i,'u.logbox'] <- u.logbox
      dtf.m.P[i,'u.bar'] <- u.bar
    }


    fwrite(dtf.m.P,file=paste0(path.main,'OUT/',name.prcp.out),sep=',')
  }
  
  if(1) # temperature residuals
  {
    setwd(paste0(path.main,'IN_temperature_residuals/'))
    list.stations <- list.files(pattern = "\\.csv$")
    n.station <- length(list.stations)
    dtf.m.T <- data.frame(stations=substr(list.stations,1,11),n=rep(NA,n.station),kim=rep(NA,n.station),hub=rep(NA,n.station),sch=rep(NA,n.station),ley=rep(NA,n.station),logbox=rep(NA,n.station),bar=rep(NA,n.station),l.kim=rep(NA,n.station),l.hub=rep(NA,n.station),l.sch=rep(NA,n.station),l.ley=rep(NA,n.station),l.logbox=rep(NA,n.station),l.bar=rep(NA,n.station),u.kim=rep(NA,n.station),u.hub=rep(NA,n.station),u.sch=rep(NA,n.station),u.ley=rep(NA,n.station),u.logbox=rep(NA,n.station),u.bar=rep(NA,n.station),stringsAsFactors = F)
    for(i in 1:n.station)
    {
      print(c(i,n.station))
      d <- read.table(list.stations[i],sep=',',header=T,stringsAsFactors = F)
      y5 <- d[,1]
      y10 <- d[,2]
      y20 <- d[,3]
      y5 <- y5[!is.na(y5)]
      y10 <- y10[!is.na(y10)]
      y20 <- y20[!is.na(y20)]
      
      n5 <- length(y5)
      n10 <- length(y10)
      n20 <- length(y20)
    
      
      n.all <- round(mean(c(n5,n10,n20)),digits=4)
      n.out.kim <- c()
      n.out.hub <- c()
      n.out.sch <- c()
      n.out.ley <- c()
      n.out.logbox <- c()
      n.out.bar <- c()
      
      l.kim <- c()
      l.hub <- c()
      l.sch <- c()
      l.ley <- c()
      l.logbox <- c()
      l.bar <- c()
      
      u.kim <- c()
      u.hub <- c()
      u.sch <- c()
      u.ley <- c()
      u.logbox <- c()
      u.bar <- c()
      
      for(j in 1:3)
      {
        if(j==1)
        {
          d.loop <- y5
          n.loop <- n5
        }
        if(j==2)
        {
          d.loop <- y10
          n.loop <- n10
        }
        if(j==3)
        {
          d.loop <- y20
          n.loop <- n20
        }
        
        ###
        MC <- as.numeric(medcouple(d.loop))
        MAD0 <- as.numeric(mad(d.loop))
        q0.125 <- as.numeric(quantile(d.loop,0.125))
        q0.25 <- as.numeric(quantile(d.loop,0.25))
        q0.375 <- as.numeric(quantile(d.loop,0.375))
        q0.50 <- as.numeric(quantile(d.loop,0.50))
        q0.625 <- as.numeric(quantile(d.loop,0.625))
        q0.75 <- as.numeric(quantile(d.loop,0.75))
        q0.875 <- as.numeric(quantile(d.loop,0.875))
        coeff.outlier <- 'auto'
        
        nlu.kim <- f.Kimber(d.loop,q0.25,q0.50,q0.75)
        nlu.hub <- f.Hubert(d.loop,q0.25,q0.75,MC)
        nlu.sch <- f.Schwertman(d.loop,q0.25,q0.50,q0.75,Z.schwertman,k.n.Schwertman[11])
        nlu.ley <- f.Leys(d.loop,q0.50,MAD0)
        nlu.logbox <- f.LogBox(d.loop,coeff.outlier,q0.125,q0.25,q0.375,q0.625,q0.75,q0.875)
        nlu.bar <- f.Barbato(d.loop,q0.25,q0.75)
        
        
        n.out.kim <- c(n.out.kim,nlu.kim[1])
        n.out.hub <- c(n.out.hub,nlu.hub[1])
        n.out.sch <- c(n.out.sch,nlu.sch[1])
        n.out.ley <- c(n.out.ley,nlu.ley[1])
        n.out.logbox <- c(n.out.logbox,nlu.logbox[1])
        n.out.bar <- c(n.out.bar,nlu.bar[1])
        
        l.kim <- c(l.kim,nlu.kim[2])
        l.hub <- c(l.hub,nlu.hub[2])
        l.sch <- c(l.sch,nlu.sch[2])
        l.ley <- c(l.ley,nlu.ley[2])
        l.logbox <- c(l.logbox,nlu.logbox[2])
        l.bar <- c(l.bar,nlu.bar[2])
        
        u.kim <- c(u.kim,nlu.kim[3])
        u.hub <- c(u.hub,nlu.hub[3])
        u.sch <- c(u.sch,nlu.sch[3])
        u.ley <- c(u.ley,nlu.ley[3])
        u.logbox <- c(u.logbox,nlu.logbox[3])
        u.bar <- c(u.bar,nlu.bar[3])
      }
      
      n.out.kim <- round(mean(n.out.kim),digits=4)
      n.out.hub <- round(mean(n.out.hub),digits=4)
      n.out.sch <- round(mean(n.out.sch),digits=4)
      n.out.ley <- round(mean(n.out.ley),digits=4)
      n.out.logbox <- round(mean(n.out.logbox),digits=4)
      n.out.bar <- round(mean(n.out.bar),digits=4)
      
      l.kim <- round(mean(l.kim),digits=4)
      l.hub <- round(mean(l.hub),digits=4)
      l.sch <- round(mean(l.sch),digits=4)
      l.ley <- round(mean(l.ley),digits=4)
      l.logbox <- round(mean(l.logbox),digits=4)
      l.bar <- round(mean(l.bar),digits=4)
      
      u.kim <- round(mean(u.kim),digits=4)
      u.hub <- round(mean(u.hub),digits=4)
      u.sch <- round(mean(u.sch),digits=4)
      u.ley <- round(mean(u.ley),digits=4)
      u.logbox <- round(mean(u.logbox),digits=4)
      u.bar <- round(mean(u.bar),digits=4)
      
      dtf.m.T[i,'n'] <- n.all
      dtf.m.T[i,'kim'] <- n.out.kim
      dtf.m.T[i,'hub'] <- n.out.hub
      dtf.m.T[i,'sch'] <- n.out.sch
      dtf.m.T[i,'ley'] <- n.out.ley
      dtf.m.T[i,'logbox'] <- n.out.logbox
      dtf.m.T[i,'bar'] <- n.out.bar
      
      dtf.m.T[i,'l.kim'] <- l.kim
      dtf.m.T[i,'l.hub'] <- l.hub
      dtf.m.T[i,'l.sch'] <- l.sch
      dtf.m.T[i,'l.ley'] <- l.ley
      dtf.m.T[i,'l.logbox'] <- l.logbox
      dtf.m.T[i,'l.bar'] <- l.bar
      
      dtf.m.T[i,'u.kim'] <- u.kim
      dtf.m.T[i,'u.hub'] <- u.hub
      dtf.m.T[i,'u.sch'] <- u.sch
      dtf.m.T[i,'u.ley'] <- u.ley
      dtf.m.T[i,'u.logbox'] <- u.logbox
      dtf.m.T[i,'u.bar'] <- u.bar
    }
    
    fwrite(dtf.m.T,file=paste0(path.main,'OUT/',name.temp.out),sep=',')
  }
}

# create histograms for all stations
if(1)
{
  count.hist <- function(x,seq.x,dtx)
  {
    x <- x[min(seq.x) <= x & x <= max(seq.x)]
    x.counts <- as.numeric(hist(x,breaks=seq.x,plot=FALSE)$counts)

    dtx[,2] <- dtx[,2]+x.counts/3 # divided by 3 because there are 3 replicas
    return(dtx)
  }
  
  if(1)
  {
  min.hist <- -500
  max.hist <- 2000
  break.hist <- 10
  
  seq.hist <- seq(from = min.hist-break.hist/2,to=max.hist+break.hist/2,by=break.hist)
  nn <- length(seq.hist)
  data.sum.hist <- data.frame(seq.hist=(seq.hist[2:nn]+seq.hist[1:(nn-1)])/2,n.points=rep(0,nn-1))
  

  setwd(paste0(path.main,'IN_precipitation_residuals/'))
  list.stations <- list.files(pattern = "\\.csv$")
  n.station <- length(list.stations)
  for(i in 1:n.station)
  {
    print(c(i,n.station))
    d <- read.table(list.stations[i],sep=',',header=T,stringsAsFactors = F)
    y5 <- d[,1]
    y10 <- d[,2]
    y20 <- d[,3]
    y5 <- y5[!is.na(y5)]
    y10 <- y10[!is.na(y10)]
    y20 <- y20[!is.na(y20)]
    
    data.sum.hist <-  count.hist(y5,seq.hist,data.sum.hist)
    data.sum.hist <-  count.hist(y10,seq.hist,data.sum.hist)
    data.sum.hist <-  count.hist(y20,seq.hist,data.sum.hist)
  }
  
  data.sum.hist <- data.sum.hist[data.sum.hist[,2] != 0,]
  n.median.P <- median(data.sum.hist[,2],na.rm=T)
  # data.sum.hist <- data.sum.hist[data.sum.hist[,2] >= 100,]

  fwrite(data.sum.hist,file=paste0(path.main,'OUT/',name.hist.prcp.out),sep=',')
  }
  
  if(1)
  {
    min.hist <- -100
    max.hist <- 100
    break.hist <- 1
    
    seq.hist <- seq(from = min.hist-break.hist/2,to=max.hist+break.hist/2,by=break.hist)
    nn <- length(seq.hist)
    data.sum.hist <- data.frame(seq.hist=(seq.hist[2:nn]+seq.hist[1:(nn-1)])/2,n.points=rep(0,nn-1))
    
    
    setwd(paste0(path.main,'IN_temperature_residuals/'))
    list.stations <- list.files(pattern = "\\.csv$")
    n.station <- length(list.stations)
    for(i in 1:n.station)
    {
      print(c(i,n.station))
      d <- read.table(list.stations[i],sep=',',header=T,stringsAsFactors = F)
      y5 <- d[,1]
      y10 <- d[,2]
      y20 <- d[,3]
      y5 <- y5[!is.na(y5)]
      y10 <- y10[!is.na(y10)]
      y20 <- y20[!is.na(y20)]
      
      data.sum.hist <-  count.hist(y5,seq.hist,data.sum.hist)
      data.sum.hist <-  count.hist(y10,seq.hist,data.sum.hist)
      data.sum.hist <-  count.hist(y20,seq.hist,data.sum.hist)
    }
    
    data.sum.hist <- data.sum.hist[data.sum.hist[,2] != 0,]
    n.median.T <- median(data.sum.hist[,2],na.rm=T)
    # data.sum.hist <- data.sum.hist[data.sum.hist[,2] >= 100,]
    
    fwrite(data.sum.hist,file=paste0(path.main,'OUT/',name.hist.temp.out),sep=',')
  }
}

# calculate the number of points removed for each station (sample size varying from n = 10 to n = 100)
if(1)
{
  if(1) # precipitation residuals
  {
    setwd(paste0(path.main,'IN_precipitation_residuals/'))
    list.stations <- list.files(pattern = "\\.csv$")
    n.station <- length(list.stations)
    seq.samples <- seq(from=10,to=100,by=10)
    n.times <- length(seq.samples)
    n.all <- n.times*n.station
    dt.P.small <- data.frame(stations=rep('AAA',n.all),n=rep(NA,n.all),n.samples=rep(NA,n.all),kim=rep(NA,n.all),hub=rep(NA,n.all),sch=rep(NA,n.all),ley=rep(NA,n.all),logbox=rep(NA,n.all),bar=rep(NA,n.all),stringsAsFactors = F)
    
    k.loop <- 0
    for(i in 1:n.station)
    {
      print(c(i,n.station))
      d <- read.table(list.stations[i],sep=',',header=T,stringsAsFactors = F)
      y5 <- d[,1]
      y10 <- d[,2]
      y20 <- d[,3]
      y5 <- y5[!is.na(y5)]
      y10 <- y10[!is.na(y10)]
      y20 <- y20[!is.na(y20)]
      
      n5 <- length(y5)
      n10 <- length(y10)
      n20 <- length(y20)
      
      
      if(min(c(n5,n10,n20)) >= max(1000,max(seq.samples)))
      {
        for(k.samples in 1:length(seq.samples))
        {
          k.loop <- k.loop+1
          dt.P.small[k.loop,'stations'] <- substr(list.stations[i],1,11)
          dt.P.small[k.loop,'n.samples'] <- seq.samples[k.samples]
          
          n.blocks <- ceiling(1000/seq.samples[k.samples])
          n.out.kim <- c()
          n.out.hub <- c()
          n.out.sch <- c()
          n.out.ley <- c()
          n.out.logbox <- c()
          n.out.bar <- c()
          
          for(j in 1:3)
          {
            if(j==1)
            {
              d.loop <- y5
            }
            if(j==2)
            {
              d.loop <- y10
            }
            if(j==3)
            {
              d.loop <- y20
            }
            
            for(l in 1:n.blocks)
            {
              d.samp <- sample(d.loop,seq.samples[k.samples])
              
              ###
              MC <- as.numeric(medcouple(d.samp))
              MAD0 <- as.numeric(mad(d.samp))
              q0.125 <- as.numeric(quantile(d.samp,0.125))
              q0.25 <- as.numeric(quantile(d.samp,0.25))
              q0.375 <- as.numeric(quantile(d.samp,0.375))
              q0.50 <- as.numeric(quantile(d.samp,0.50))
              q0.625 <- as.numeric(quantile(d.samp,0.625))
              q0.75 <- as.numeric(quantile(d.samp,0.75))
              q0.875 <- as.numeric(quantile(d.samp,0.875))
              coeff.outlier <- 'auto'
              
              nlu.kim <- f.Kimber(d.samp,q0.25,q0.50,q0.75)
              nlu.hub <- f.Hubert(d.samp,q0.25,q0.75,MC)
              nlu.sch <- f.Schwertman(d.samp,q0.25,q0.50,q0.75,Z.schwertman,k.n.Schwertman[11])
              nlu.ley <- f.Leys(d.samp,q0.50,MAD0)
              nlu.logbox <- f.LogBox(d.samp,coeff.outlier,q0.125,q0.25,q0.375,q0.625,q0.75,q0.875)
              nlu.bar <- f.Barbato(d.samp,q0.25,q0.75)
              
              
              n.out.kim <- c(n.out.kim,nlu.kim[1])
              n.out.hub <- c(n.out.hub,nlu.hub[1])
              n.out.sch <- c(n.out.sch,nlu.sch[1])
              n.out.ley <- c(n.out.ley,nlu.ley[1])
              n.out.logbox <- c(n.out.logbox,nlu.logbox[1])
              n.out.bar <- c(n.out.bar,nlu.bar[1])
              
            }
          }
          
          
          n.out.kim <- round(sum(n.out.kim)/3,digits=4)
          n.out.hub <- round(sum(n.out.hub)/3,digits=4)
          n.out.sch <- round(sum(n.out.sch)/3,digits=4)
          n.out.ley <- round(sum(n.out.ley)/3,digits=4)
          n.out.logbox <- round(sum(n.out.logbox)/3,digits=4)
          n.out.bar <- round(sum(n.out.bar)/3,digits=4)
          
          dt.P.small[k.loop,'n'] <- seq.samples[k.samples]*n.blocks
          dt.P.small[k.loop,'kim'] <- n.out.kim
          dt.P.small[k.loop,'hub'] <- n.out.hub
          dt.P.small[k.loop,'sch'] <- n.out.sch
          dt.P.small[k.loop,'ley'] <- n.out.ley
          dt.P.small[k.loop,'logbox'] <- n.out.logbox
          dt.P.small[k.loop,'bar'] <- n.out.bar
          
          
        }
      }
    }
    
    dt.P.small <- dt.P.small[!is.na(dt.P.small[,'n']),]

    fwrite(dt.P.small,file=paste0(path.main,'OUT/',name.small.prcp.out),sep=',')
  }
  
  if(1) # temperature residuals
  {
    setwd(paste0(path.main,'IN_temperature_residuals/'))
    list.stations <- list.files(pattern = "\\.csv$")
    n.station <- length(list.stations)
    seq.samples <- seq(from=10,to=100,by=10)
    n.times <- length(seq.samples)
    n.all <- n.times*n.station
    dt.T.small <- data.frame(stations=rep('AAA',n.all),n=rep(NA,n.all),n.samples=rep(NA,n.all),kim=rep(NA,n.all),hub=rep(NA,n.all),sch=rep(NA,n.all),ley=rep(NA,n.all),logbox=rep(NA,n.all),bar=rep(NA,n.all),stringsAsFactors = F)
    
    k.loop <- 0
    for(i in 1:n.station)
    {
      print(c(i,n.station))
      d <- read.table(list.stations[i],sep=',',header=T,stringsAsFactors = F)
      y5 <- d[,1]
      y10 <- d[,2]
      y20 <- d[,3]
      y5 <- y5[!is.na(y5)]
      y10 <- y10[!is.na(y10)]
      y20 <- y20[!is.na(y20)]
      
      n5 <- length(y5)
      n10 <- length(y10)
      n20 <- length(y20)
      
      
      if(min(c(n5,n10,n20)) >= 1000)
      {
        for(k.samples in 1:length(seq.samples))
        {
          k.loop <- k.loop+1
          dt.T.small[k.loop,'stations'] <- substr(list.stations[i],1,11)
          dt.T.small[k.loop,'n.samples'] <- seq.samples[k.samples]
          
          n.blocks <- ceiling(1000/seq.samples[k.samples])
          n.out.kim <- c()
          n.out.hub <- c()
          n.out.sch <- c()
          n.out.ley <- c()
          n.out.logbox <- c()
          n.out.bar <- c()
          
          for(j in 1:3)
          {
            if(j==1)
            {
              d.loop <- y5
            }
            if(j==2)
            {
              d.loop <- y10
            }
            if(j==3)
            {
              d.loop <- y20
            }
            
            for(l in 1:n.blocks)
            {
              d.samp <- sample(d.loop,seq.samples[k.samples])
              
              ###
              MC <- as.numeric(medcouple(d.samp))
              MAD0 <- as.numeric(mad(d.samp))
              q0.125 <- as.numeric(quantile(d.samp,0.125))
              q0.25 <- as.numeric(quantile(d.samp,0.25))
              q0.375 <- as.numeric(quantile(d.samp,0.375))
              q0.50 <- as.numeric(quantile(d.samp,0.50))
              q0.625 <- as.numeric(quantile(d.samp,0.625))
              q0.75 <- as.numeric(quantile(d.samp,0.75))
              q0.875 <- as.numeric(quantile(d.samp,0.875))
              coeff.outlier <- 'auto'
              
              nlu.kim <- f.Kimber(d.samp,q0.25,q0.50,q0.75)
              nlu.hub <- f.Hubert(d.samp,q0.25,q0.75,MC)
              nlu.sch <- f.Schwertman(d.samp,q0.25,q0.50,q0.75,Z.schwertman,k.n.Schwertman[11])
              nlu.ley <- f.Leys(d.samp,q0.50,MAD0)
              nlu.logbox <- f.LogBox(d.samp,coeff.outlier,q0.125,q0.25,q0.375,q0.625,q0.75,q0.875)
              nlu.bar <- f.Barbato(d.samp,q0.25,q0.75)
              
              n.out.kim <- c(n.out.kim,nlu.kim[1])
              n.out.hub <- c(n.out.hub,nlu.hub[1])
              n.out.sch <- c(n.out.sch,nlu.sch[1])
              n.out.ley <- c(n.out.ley,nlu.ley[1])
              n.out.logbox <- c(n.out.logbox,nlu.logbox[1])
              n.out.bar <- c(n.out.bar,nlu.bar[1])
            }
          }
          
          
          n.out.kim <- round(sum(n.out.kim)/3,digits=4)
          n.out.hub <- round(sum(n.out.hub)/3,digits=4)
          n.out.sch <- round(sum(n.out.sch)/3,digits=4)
          n.out.ley <- round(sum(n.out.ley)/3,digits=4)
          n.out.logbox <- round(sum(n.out.logbox)/3,digits=4)
          n.out.bar <- round(sum(n.out.bar)/3,digits=4)
          
          dt.T.small[k.loop,'n'] <- seq.samples[k.samples]*n.blocks
          dt.T.small[k.loop,'kim'] <- n.out.kim
          dt.T.small[k.loop,'hub'] <- n.out.hub
          dt.T.small[k.loop,'sch'] <- n.out.sch
          dt.T.small[k.loop,'ley'] <- n.out.ley
          dt.T.small[k.loop,'logbox'] <- n.out.logbox
          dt.T.small[k.loop,'bar'] <- n.out.bar
          
          
        }
      }
    }
    
    dt.T.small <- dt.T.small[!is.na(dt.T.small[,'n']),]
    
    fwrite(dt.T.small,file=paste0(path.main,'OUT/',name.small.temp.out),sep=',')
  }
}











