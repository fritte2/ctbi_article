setwd("~/Science/ctbi3/hess-2021-609-main_v1.0.0/PartII.aggregation/step2_contamination.of.datasets")

set.seed(1)
library(ctbi)
library(forecast)
# these two functions randomly create missing values and gap in the raw data.
source('create.poisson.blocks_v1.R')
source('create.3big.blocks_v1.R')

data.temperature <- read.table('temperature_SJER.csv',header=T,sep=',',stringsAsFactors = F)
data.precipitation <- read.table('PRCP_CapeLeuwin.csv',header=T,sep=',',stringsAsFactors = F)
data.CH4 <- read.table('CH4_epica.csv',header=T,sep=',',stringsAsFactors = F)

data.temperature[,1] <- as.POSIXct(data.temperature[,1],tz='UTC')
data.precipitation[,1] <- as.Date(data.precipitation[,1])

# replace the rare missing values in the raw data with interpolated values
interp.temp <- approx(x=as.numeric(data.temperature[,1]),y=data.temperature[,2],xout=as.numeric(data.temperature[,1]))
data.temperature[,2] <- interp.temp[[2]]
interp.temp <- approx(x=as.numeric(data.precipitation[,1]),y=data.precipitation[,2],xout=as.numeric(data.precipitation[,1]))
data.precipitation[,2] <- interp.temp[[2]]
interp.temp <- approx(x=as.numeric(data.CH4[,1]),y=data.CH4[,2],xout=as.numeric(data.CH4[,1]))
data.CH4[,2] <- interp.temp[[2]]

# cut for a better visualization
read1 <- as.POSIXct('2020-08-01') <= data.temperature[,1]
read2 <- as.Date('1990-01-01') <= data.precipitation[,1]
data.temperature <- data.temperature[read1,]
data.precipitation <- data.precipitation[read2,]

# input for all the datasets:
# contamination inputs : 30% of the dataset is contaminated
frac.NA.blocks <- 0.2 # 20% of large data gaps
frac.NA <- 0.095 # 9.5% of missing values
frac.outliers <- 0.005 # 0.5% of outliers

# aggregate the raw and contaminated datasets.
k.outliers <- 0.6
for(k.datasets in 1:3)
{
  if(k.datasets==1)
  {
    data <- data.temperature
    bin.period <- '1 hour'
    bin.FUN <- 'mean'
    bin.center <- as.POSIXct('2020-07-01 00:00:00',tz='UTC')
    bin.side <- NULL
    SCI.min <- Inf
    ylim <- c(-Inf,Inf)
    bin.max.f.NA <- 0.2
  }
  if(k.datasets==2)
  {
    data <- data.precipitation
    bin.period <- '1 month'
    bin.FUN <- 'sum'
    bin.center <- NULL
    bin.side <- as.Date('1907-01-01')
    SCI.min <- Inf
    ylim <- c(0,Inf)
    bin.max.f.NA <- 0.2
  }
  if(k.datasets==3)
  {
    data <- data.CH4
    bin.period <- 2000
    bin.FUN <- 'mean'
    bin.center <- NULL
    bin.side <- 0
    SCI.min <- Inf
    ylim <- c(-Inf,Inf)
    bin.max.f.NA <- 1
  }

  N <- length(data[,1])
  data.raw <- data
  data <- data.frame(data,index.data=1:N,stringsAsFactors = F)

  N.NA.blocks <- round(frac.NA.blocks*N,digits=0)
  N.NA <- round(frac.NA*N,digits=0)
  N.outliers <- round(frac.outliers*N,digits=0)

  mean0 <- mean(data[,2],na.rm=T)
  max0 <- max(data[,2],na.rm=T)
  min0 <- min(data[,2],na.rm=T)

  index.blocks <- create.3big.blocks(data,N.NA.blocks)
  data <- data[-index.blocks,]

  index.NA <- create.poisson.blocks(data,N,N.NA)
  data <- data[-index.NA,]

  index.outliers <- create.poisson.blocks(data,N,N.outliers)
  index.min <- sample(index.outliers,round(floor(length(index.outliers)/2)))
  index.max <- setdiff(index.outliers,index.min)

  # put the outliers
  data <- data.raw
  data[index.NA,2] <- NA
  if(k.datasets != 2)
  {
    data[index.min,2] <- min0-abs(mean0-min0)*0.5
    data[index.max,2] <- max0+abs(mean0-max0)*0.5
  }else
  {
    data[index.min,2] <- max0+abs(mean0-max0)*0.5
    data[index.max,2] <- max0+abs(mean0-max0)*0.5
  }

  ylim.min <- c(5,0,200)
  ylim.max <- c(50,300,1200)
  index.out <- c(index.min,index.max)
  # plot of the data #
  if(k.datasets != 3)
  {
    d.t <- data
    d.t[index.blocks,2] <- NA # time gap

    if(k.datasets==1)
    {
      plot(d.t[,1],d.t[,2],ylim=c(ylim.min[k.datasets],ylim.max[k.datasets]),xlim=c(min(data.raw[,1]),max(data.raw[,1])),type='l')
    }
    if(k.datasets==2)
    {
      plot(d.t[,1],d.t[,2],ylim=c(ylim.min[k.datasets],ylim.max[k.datasets]),xlim=c(min(data.raw[,1]),max(data.raw[,1])),type='h')
    }
  }else
  {
    d.t <- data
    d.t[,1] <- -d.t[,1]
    d.t[index.blocks,2] <- NA

    plot(d.t[,1],d.t[,2],ylim=c(ylim.min[k.datasets],ylim.max[k.datasets]),xlim=c(-max(data.raw[,1]),-min(data.raw[,1])),type='l')
  }

  data <- data[-index.blocks,]


  # plot of the raw data
  if(k.datasets != 3)
  {
    if(k.datasets==1)
    {
      plot(data.raw[,1],data.raw[,2],type='l',ylim=c(ylim.min[k.datasets],ylim.max[k.datasets]))
    }
    if(k.datasets==2)
    {
      plot(data.raw[,1],data.raw[,2],type='h',ylim=c(ylim.min[k.datasets],ylim.max[k.datasets]))
    }
  }else
  {
    d.t <- data.raw
    d.t[,1] <- -d.t[,1]
    plot(d.t[,1],d.t[,2],type='l',ylim=c(ylim.min[k.datasets],ylim.max[k.datasets]))
  }


  list.main.raw <- ctbi(data.raw,bin.side=bin.side,bin.period=bin.period,bin.center=bin.center,bin.FUN = bin.FUN,bin.max.f.NA = bin.max.f.NA,SCI.min = SCI.min,ylim=ylim,k.outliers=Inf)
  list.main <- ctbi(data,bin.side=bin.side,bin.period=bin.period,bin.center=bin.center,bin.FUN = bin.FUN,bin.max.f.NA = bin.max.f.NA,SCI.min = SCI.min,ylim=ylim,k.outliers=k.outliers)

  data0.raw <- list.main.raw$data0
  data0 <- list.main$data0

  data1.raw <- list.main.raw$data1
  data1 <- list.main$data1

  # save the original data.set
  if(k.datasets==1)
  {
    temperature.aggreg0 <- data0
    temperature.aggreg0.raw <- data0.raw
    temperature.aggreg1 <- data1
    temperature.aggreg1.raw <- data1.raw
  }
  if(k.datasets==2)
  {
    precipitation.aggreg0 <- data0
    precipitation.aggreg0.raw <- data0.raw
    precipitation.aggreg1 <- data1
    precipitation.aggreg1.raw <- data1.raw
  }
  if(k.datasets==3)
  {
    CH4.aggreg0 <- data0
    CH4.aggreg0.raw <- data0.raw
    CH4.aggreg1 <- data1
    CH4.aggreg1.raw <- data1.raw
  }


  # second aggregation
  if(k.datasets==1)
  {
    data <- data1[,c(1,2)]
    data.raw <- data1.raw[,c(1,2)]
    bin.period <- '1 day'
    bin.FUN <- 'mean'
    bin.center <- NULL
    bin.side <- as.POSIXct('2020-07-01 00:00:00',tz='UTC')
    SCI.min <- 0.6
    ylim <- c(-Inf,Inf)
    bin.max.f.NA <- 0.2
  }
  if(k.datasets==2)
  {
    data <- data1[,c(1,2)]
    data.raw <- data1.raw[,c(1,2)]
    bin.period <- '1 year'
    bin.FUN <- 'sum'
    bin.center <- NULL
    bin.side <- as.Date('1907-01-01')
    SCI.min <- 0.6 # inf if no imputation
    ylim <- c(0,Inf)
    bin.max.f.NA <- 0.2
  }
  if(k.datasets==3)
  {
    data <- data1[,c(1,2)]
    data.raw <- data1.raw[,c(1,2)]
    bin.period <- 20000
    bin.FUN <- 'mean'
    bin.center <- NULL
    bin.side <- 0
    SCI.min <- 0.6
    ylim <- c(-Inf,Inf)
    bin.max.f.NA <- 1
  }

  list.main.raw <- ctbi(data.raw,bin.side=bin.side,bin.period=bin.period,bin.center=bin.center,bin.FUN = bin.FUN,bin.max.f.NA = bin.max.f.NA,SCI.min = SCI.min,ylim=ylim,k.outliers=Inf)

  list.main <- ctbi(data,bin.side=bin.side,bin.period=bin.period,bin.center=bin.center,bin.FUN = bin.FUN,bin.max.f.NA = bin.max.f.NA,SCI.min = SCI.min,ylim=ylim,k.outliers=k.outliers)

  data1.raw <- list.main.raw$data1
  data1 <- list.main$data1

  SCI.raw <- list.main.raw$SCI
  SCI <- list.main$SCI



  mean.cycle.raw <- list.main.raw$mean.cycle
  mean.cycle <- list.main$mean.cycle

  # save the original data.set
  if(k.datasets==1)
  {
    temperature.aggreg2 <- data1
    temperature.aggreg2.raw <- data1.raw

    mean.cycle.temperature <- mean.cycle
    mean.cycle.temperature.raw <- mean.cycle.raw

    SCI.temperature <- SCI
    SCI.temperature.raw <- SCI.raw
  }
  if(k.datasets==2)
  {
    precipitation.aggreg2 <- data1
    precipitation.aggreg2.raw <- data1.raw

    mean.cycle.precipitation <- mean.cycle
    mean.cycle.precipitation.raw <- mean.cycle.raw

    SCI.precipitation <- SCI
    SCI.precipitation.raw <- SCI.raw
  }
  if(k.datasets==3)
  {
    CH4.aggreg2 <- data1
    CH4.aggreg2.raw <- data1.raw

    mean.cycle.CH4 <- mean.cycle
    mean.cycle.CH4.raw <- mean.cycle.raw

    SCI.CH4 <- SCI
    SCI.CH4.raw <- SCI.raw
  }
}


# calculate the false positive and false negative for this study
y.outlier.min <- c(10,-1,255)
y.outlier.max <- c(47,255,1050)
for(k.datasets in 1:3)
{
  if(k.datasets==1)
  {
    data0 <- temperature.aggreg0
    data0.raw <- temperature.aggreg0.raw
    data1 <- temperature.aggreg1
    data1.raw <- temperature.aggreg1.raw
  }
  if(k.datasets==2)
  {
    data0 <- precipitation.aggreg0
    data0.raw <- precipitation.aggreg0.raw
    data1 <- precipitation.aggreg1
    data1.raw <- precipitation.aggreg1.raw
  }
  if(k.datasets==3)
  {
    data0 <- CH4.aggreg0
    data0.raw <- CH4.aggreg0.raw
    data1 <- CH4.aggreg1
    data1.raw <- CH4.aggreg1.raw
  }

  read.good.bin <- (data0[,'index.bin'] > 0)

  data0 <- data0[read.good.bin,]
  data0.raw <- data0.raw[read.good.bin,]

  N0 <- length(data0[,1])
  index0 <- 1:N0

  d.temp <- data0[,2]
  d.temp[!is.na(data0[,'outliers'])] <- data0[!is.na(data0[,'outliers']),'outliers']
  
  print(length(d.temp))
  print(N0)

  read.true <- d.temp < y.outlier.min[k.datasets] | y.outlier.max[k.datasets] < d.temp
  read.true[is.na(read.true)] <- F

  index.outliers <- index0[read.true]
  index.flagged <- index0[!is.na(data0[,'outliers'])]

  index.true.positive <- intersect(index.flagged,index.outliers)

  N.true.positive <- length(index.true.positive)
  N.false.positive <- length(index.flagged)-N.true.positive

  index.false.negative <- setdiff(index.outliers,index.flagged)
  N.false.negative <- length(index.false.negative)

  index.not.outliers <- setdiff(index0,index.outliers)
  index.not.flagged <- index0[is.na(data0[,'outliers'])]

  index.true.negative <- intersect(index.not.flagged,index.not.outliers)
  N.true.negative <- length(index.true.negative)

  print('#####')
  if(k.datasets == 1)
  {
    print('temperature dataset:')
  }
  if(k.datasets == 2)
  {
    print('precipitation dataset:')
  }
  if(k.datasets == 3)
  {
    print('CH4 dataset:')
  }
  # print(sum(data1.raw[,'index.bin'] < 0))
  # print(sum(data1[,'index.bin'] < 0))

  print(paste0('true positive: ',N.true.positive))
  print(paste0('false positive: ',N.false.positive))
  print(paste0('true negative: ',N.true.negative))
  print(paste0('false negative: ',N.false.negative))
}


# calculate the false positive and false negative with tsoutliers
y.outlier.min <- c(10,-1,255)
y.outlier.max <- c(47,255,1050)
for(k.datasets in 1:3)
{
  if(k.datasets==1)
  {
    data0 <- temperature.aggreg0
    data0.raw <- temperature.aggreg0.raw
  }
  if(k.datasets==2)
  {
    data0 <- precipitation.aggreg0
    data0.raw <- precipitation.aggreg0.raw
  }
  if(k.datasets==3)
  {
    data0 <- CH4.aggreg0
    data0.raw <- CH4.aggreg0.raw
  }
  
  N0 <- length(data0[,1])
  index0 <- 1:N0
  
  d.temp <- data0[,2]
  d.temp[!is.na(data0[,'outliers'])] <- data0[!is.na(data0[,'outliers']),'outliers'] # replace flagged outlier in the timeseries
  
  read.true <- d.temp < y.outlier.min[k.datasets] | y.outlier.max[k.datasets] < d.temp
  read.true[is.na(read.true)] <- F
  
  index.outliers <- index0[read.true]
  # index.flagged <- tsoutliers(as.numeric(d.temp),lambda='auto')$index
  index.flagged <- tsoutliers(as.numeric(d.temp),lambda=NULL)$index
  
  index.true.positive <- intersect(index.flagged,index.outliers)
  
  N.true.positive <- length(index.true.positive)
  N.false.positive <- length(index.flagged)-N.true.positive
  
  index.false.negative <- setdiff(index.outliers,index.flagged)
  N.false.negative <- length(index.false.negative)
  
  index.not.outliers <- setdiff(index0,index.outliers)
  index.not.flagged <- setdiff(index0,index.flagged)
  
  index.true.negative <- intersect(index.not.flagged,index.not.outliers)
  N.true.negative <- length(index.true.negative)
  
  print('#####')
  if(k.datasets == 1)
  {
    print('tsoutliers, temperature dataset:')
  }
  if(k.datasets == 2)
  {
    print('tsoutliers, precipitation dataset:')
  }
  if(k.datasets == 3)
  {
    print('tsoutliers, CH4 dataset:')
  }
  # print(sum(data1.raw[,'index.bin'] < 0))
  # print(sum(data1[,'index.bin'] < 0))
  
  print(paste0('true positive: ',N.true.positive))
  print(paste0('false positive: ',N.false.positive))
  print(paste0('true negative: ',N.true.negative))
  print(paste0('false negative: ',N.false.negative))
}


# plotting
ylim.min <- c(15,0,350)
ylim.max <- c(40,300,800)
plot(temperature.aggreg1.raw[,1],temperature.aggreg1.raw[,2],type='l',ylim=c(ylim.min[1],ylim.max[1]))
lines(temperature.aggreg1[,1],temperature.aggreg1[,2],col='red')
lines(temperature.aggreg2.raw[,1],temperature.aggreg2.raw[,2],lwd=3)
lines(temperature.aggreg2[,1],temperature.aggreg2[,2],col='red',lwd=3)
plot(precipitation.aggreg1.raw[,1],precipitation.aggreg1.raw[,2],type='l',ylim=c(ylim.min[2],ylim.max[2]))
lines(precipitation.aggreg1[,1],precipitation.aggreg1[,2],col='red')

CH4.aggreg1.raw[,1] <- -CH4.aggreg1.raw[,1]
CH4.aggreg1[,1] <- -CH4.aggreg1[,1]
CH4.aggreg2.raw[,1] <- -CH4.aggreg2.raw[,1]
CH4.aggreg2[,1] <- -CH4.aggreg2[,1]
plot(CH4.aggreg1.raw[,1],CH4.aggreg1.raw[,2],type='l',ylim=c(ylim.min[3],ylim.max[3]))
lines(CH4.aggreg1[,1],CH4.aggreg1[,2],col='red')
lines(CH4.aggreg2.raw[,1],CH4.aggreg2.raw[,2],lwd=3)
lines(CH4.aggreg2[,1],CH4.aggreg2[,2],col='red',lwd=3)

plot(mean.cycle.temperature.raw[,1],mean.cycle.temperature.raw[,2],type='l',lwd=2,ylim=c(-6,6),xlim=c(temperature.aggreg2.raw[1,3],temperature.aggreg2.raw[1,4]))
lines(mean.cycle.temperature[,1],mean.cycle.temperature[,2],col='red',lwd=2)
lines(mean.cycle.temperature[,1],mean.cycle.temperature[,2]+mean.cycle.temperature[,3],col='red',lty=2)
lines(mean.cycle.temperature[,1],mean.cycle.temperature[,2]-mean.cycle.temperature[,3],col='red',lty=2)

# Months are centered around day ~16. Put the day to 1 for a better visualization
mean.cycle.precipitation.raw[,1] <- as.Date(paste0(substr(as.character(mean.cycle.precipitation.raw[,1]),1,7),'-01'))
plot(mean.cycle.precipitation.raw[,1],mean.cycle.precipitation.raw[,2],type='l',lwd=2,ylim=c(-60,120),xlim=c(precipitation.aggreg2.raw[1,3],precipitation.aggreg2.raw[1,4]))
lines(mean.cycle.precipitation.raw[,1],mean.cycle.precipitation[,2],col='red',lwd=2)
lines(mean.cycle.precipitation.raw[,1],mean.cycle.precipitation[,2]+mean.cycle.precipitation[,3],col='red',lty=2)
lines(mean.cycle.precipitation.raw[,1],mean.cycle.precipitation[,2]-mean.cycle.precipitation[,3],col='red',lty=2)

plot(mean.cycle.CH4.raw[,1],mean.cycle.CH4.raw[,2],type='l',lwd=2,ylim=c(-80,80),xlim=c(CH4.aggreg2.raw[1,3],CH4.aggreg2.raw[1,4]))
lines(mean.cycle.CH4[,1],mean.cycle.CH4[,2],col='red',lwd=2)
lines(mean.cycle.CH4[,1],mean.cycle.CH4[,2]+mean.cycle.CH4[,3],col='red',lty=2)
lines(mean.cycle.CH4[,1],mean.cycle.CH4[,2]-mean.cycle.CH4[,3],col='red',lty=2)

per.diff.T <- mean((temperature.aggreg1[,2]-temperature.aggreg1.raw[,2])*100/temperature.aggreg1.raw[,2],na.rm=T)
per.diff.P <- mean((precipitation.aggreg1[,2]-precipitation.aggreg1.raw[,2])*100/precipitation.aggreg1.raw[,2],na.rm=T)
per.diff.CH4 <- mean((CH4.aggreg1[,2]-CH4.aggreg1.raw[,2])*100/CH4.aggreg1.raw[,2],na.rm=T)

sd.diff.T <- sd((temperature.aggreg1[,2]-temperature.aggreg1.raw[,2])*100/temperature.aggreg1.raw[,2],na.rm=T)
sd.diff.P <- sd((precipitation.aggreg1[,2]-precipitation.aggreg1.raw[,2])*100/precipitation.aggreg1.raw[,2],na.rm=T)
sd.diff.CH4 <- sd((CH4.aggreg1[,2]-CH4.aggreg1.raw[,2])*100/CH4.aggreg1.raw[,2],na.rm=T)

print('#####')
print(paste0('Temperature dataset, difference with raw data (lvl 1 of aggregation) : ',round(per.diff.T,digits=1),' +/- ',round(sd.diff.T,digits=1),' %'))
print(paste0('Precipitation dataset, difference with raw data (lvl 1 of aggregation) : ',round(per.diff.P,digits=0),' +/- ',round(sd.diff.P,digits=0),' %'))
print(paste0('CH4 dataset, difference with raw data (lvl 1 of aggregation) : ',round(per.diff.CH4,digits=1),' +/- ',round(sd.diff.CH4,digits=0),' %'))





