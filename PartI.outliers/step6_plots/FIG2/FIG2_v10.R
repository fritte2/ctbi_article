
path.main <- "~/Science/ctbi12/PartI.outliers/step6_plots/FIG2/"
setwd(path.main)

set.seed(1)

colley <- 'red'
colhub <- 'blue'
colkim <- 'green'
colsch <- 'orange'


if(1) # histograms
{
  if(1) # precipitation
  {
data.sum.hist <- read.table(paste0(path.main,'hist_precipitation_v4.csv'),sep=',',stringsAsFactors = F,header=T)

data.sum.hist <- data.sum.hist[data.sum.hist[,2] > 100,]
plot(data.sum.hist[,1],data.sum.hist[,2],type='h',lwd=4,col='black')

    
data.threshold <- read.table(paste0(path.main,'thresholds_precipitation_v4.csv'),sep=',',stringsAsFactors = F,header=T)

l.logbox <- median(data.threshold[,'l.logbox'])
u.logbox <- median(data.threshold[,'u.logbox'])

l.ley <- median(data.threshold[,'l.ley'])
u.ley <- median(data.threshold[,'u.ley'])

l.sch <- median(data.threshold[,'l.sch'])
u.sch <- median(data.threshold[,'u.sch'])

l.kim <- median(data.threshold[,'l.kim'])
u.kim <- median(data.threshold[,'u.kim'])

l.hub <- median(data.threshold[,'l.hub'])
u.hub <- median(data.threshold[,'u.hub'])


abline(v=c(l.logbox,u.logbox),lwd=2)
abline(v=c(l.ley,u.ley),col=colley,lwd=2)
abline(v=c(l.sch,u.sch),col=colsch,lwd=2)
abline(v=c(l.kim,u.kim),col=colkim,lwd=2)
abline(v=c(l.hub,u.hub),col=colhub,lwd=2)
}
  
  if(1) # temperature
  {
  data.sum.hist <- read.table(paste0(path.main,'hist_temperature_v4.csv'),sep=',',stringsAsFactors = F,header=T)
  
  data.sum.hist <- data.sum.hist[data.sum.hist[,2] > 100,]
  plot(data.sum.hist[,1],data.sum.hist[,2],type='h',lwd=4,col='black',xlim=c(-24,24))
  
  
  data.threshold <- read.table(paste0(path.main,'thresholds_temperature_v4.csv'),sep=',',stringsAsFactors = F,header=T)
  
  l.logbox <- median(data.threshold[,'l.logbox'])
  u.logbox <- median(data.threshold[,'u.logbox'])
  
  l.ley <- median(data.threshold[,'l.ley'])
  u.ley <- median(data.threshold[,'u.ley'])
  
  l.sch <- median(data.threshold[,'l.sch'])
  u.sch <- median(data.threshold[,'u.sch'])
  
  l.kim <- median(data.threshold[,'l.kim'])
  u.kim <- median(data.threshold[,'u.kim'])
  
  l.hub <- median(data.threshold[,'l.hub'])
  u.hub <- median(data.threshold[,'u.hub'])
  
  
  abline(v=c(l.logbox,u.logbox),lwd=2)
  abline(v=c(l.ley,u.ley),col=colley,lwd=2)
  abline(v=c(l.sch,u.sch),col=colsch,lwd=2)
  abline(v=c(l.kim,u.kim),col=colkim,lwd=2)
  abline(v=c(l.hub,u.hub),col=colhub,lwd=2)
  }
}


if(1) # plots for methods kim, hub, ley and sch
{
  if(1) # precipitation
  {
    dt.P <- read.table(paste0(path.main,'thresholds_precipitation_v4.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.P[,1])
    
    sd.kim <- sd(dt.P[,'kim']/dt.P[,'n'])*100
    sd.hub <- sd(dt.P[,'hub']/dt.P[,'n'])*100
    sd.sch <- sd(dt.P[,'sch']/dt.P[,'n'])*100
    sd.ley <- sd(dt.P[,'ley']/dt.P[,'n'])*100
    
    mean.kim <- mean(dt.P[,'kim']/dt.P[,'n'])*100
    sd.kim.low <- mean.kim-sd.kim
    sd.kim.up <- mean.kim+sd.kim
    
    mean.hub <- mean(dt.P[,'hub']/dt.P[,'n'])*100
    sd.hub.low <- mean.hub-sd.hub
    sd.hub.up <- mean.hub+sd.hub
    
    mean.sch <- mean(dt.P[,'sch']/dt.P[,'n'])*100
    sd.sch.low <- mean.sch-sd.sch
    sd.sch.up <- mean.sch+sd.sch
    
    mean.ley <- mean(dt.P[,'ley']/dt.P[,'n'])*100
    sd.ley.low <- mean.ley-sd.ley
    sd.ley.up <- mean.ley+sd.ley
    
    dt.P <- read.table(paste0(path.main,'thresholds_precipitation_smallsample_v6.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.P[,1])
    
    seq.sample <- seq(from=100,to=10,by=-10)
    
    x <- 120
    for(i in 1:length(seq.sample))
    {
      x <- c(x,seq.sample[i])
      read_ <- dt.P[,'n.samples'] == seq.sample[i]
      
      mean.kim.loop <- mean(dt.P[read_,'kim']/dt.P[read_,'n'])*100
      sd.kim.loop <- sd(dt.P[read_,'kim']/dt.P[read_,'n'])*100
      
      mean.kim <- c(mean.kim,mean.kim.loop)
      sd.kim.low <- c(sd.kim.low,mean.kim.loop-sd.kim.loop)
      sd.kim.up <- c(sd.kim.up,mean.kim.loop+sd.kim.loop)
      
      mean.hub.loop <- mean(dt.P[read_,'hub']/dt.P[read_,'n'])*100
      sd.hub.loop <- sd(dt.P[read_,'hub']/dt.P[read_,'n'])*100
      
      mean.hub <- c(mean.hub,mean.hub.loop)
      sd.hub.low <- c(sd.hub.low,mean.hub.loop-sd.hub.loop)
      sd.hub.up <- c(sd.hub.up,mean.hub.loop+sd.hub.loop)
      

      mean.sch.loop <- mean(dt.P[read_,'sch']/dt.P[read_,'n'])*100
      sd.sch.loop <- sd(dt.P[read_,'sch']/dt.P[read_,'n'])*100
      
      mean.sch <- c(mean.sch,mean.sch.loop)
      sd.sch.low <- c(sd.sch.low,mean.sch.loop-sd.sch.loop)
      sd.sch.up <- c(sd.sch.up,mean.sch.loop+sd.sch.loop)
      
      
      mean.ley.loop <- mean(dt.P[read_,'ley']/dt.P[read_,'n'])*100
      sd.ley.loop <- sd(dt.P[read_,'ley']/dt.P[read_,'n'])*100
      
      mean.ley <- c(mean.ley,mean.ley.loop)
      sd.ley.low <- c(sd.ley.low,mean.ley.loop-sd.ley.loop)
      sd.ley.up <- c(sd.ley.up,mean.ley.loop+sd.ley.loop)

      
    }
    
    data.out <- data.frame(sample.size=x,sd.kim.low=sd.kim.low,mean.kim=mean.kim,sd.kim.up=sd.kim.up,sd.sch.low=sd.sch.low,mean.sch=mean.sch,sd.sch.up=sd.sch.up,sd.ley.low=sd.ley.low,mean.ley=mean.ley,sd.ley.up=sd.ley.up,sd.hub.low=sd.hub.low,mean.hub=mean.hub,sd.hub.up=sd.hub.up,stringsAsFactors = F)
    data.out <- data.out[order(data.out[,'sample.size']),]
  
    
    plot(data.out[,'sample.size'],data.out[,'mean.kim'],ylim=c(3.5,12.5),col=colkim,type='l',ylab='precip')
    lines(data.out[,'sample.size'],data.out[,'sd.kim.low'],col=colkim,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.kim.up'],col=colkim,lty=2)
    
    lines(data.out[,'sample.size'],data.out[,'mean.hub'],col=colhub)
    lines(data.out[,'sample.size'],data.out[,'sd.hub.low'],col=colhub,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.hub.up'],col=colhub,lty=2)
    
    lines(data.out[,'sample.size'],data.out[,'mean.sch'],col=colsch)
    lines(data.out[,'sample.size'],data.out[,'sd.sch.low'],col=colsch,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.sch.up'],col=colsch,lty=2)
    
    lines(data.out[,'sample.size'],data.out[,'mean.ley'],col=colley)
    lines(data.out[,'sample.size'],data.out[,'sd.ley.low'],col=colley,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.ley.up'],col=colley,lty=2)
  }
  
  if(1) # temperature
  {
    dt.T <- read.table(paste0(path.main,'thresholds_temperature_v4.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.T[,1])
    
    sd.kim <- sd(dt.T[,'kim']/dt.T[,'n'])*100
    sd.hub <- sd(dt.T[,'hub']/dt.T[,'n'])*100
    sd.sch <- sd(dt.T[,'sch']/dt.T[,'n'])*100
    sd.ley <- sd(dt.T[,'ley']/dt.T[,'n'])*100
    
    mean.kim <- mean(dt.T[,'kim']/dt.T[,'n'])*100
    sd.kim.low <- mean.kim-sd.kim
    sd.kim.up <- mean.kim+sd.kim
    
    mean.hub <- mean(dt.T[,'hub']/dt.T[,'n'])*100
    sd.hub.low <- mean.hub-sd.hub
    sd.hub.up <- mean.hub+sd.hub
    
    mean.sch <- mean(dt.T[,'sch']/dt.T[,'n'])*100
    sd.sch.low <- mean.sch-sd.sch
    sd.sch.up <- mean.sch+sd.sch
    
    mean.ley <- mean(dt.T[,'ley']/dt.T[,'n'])*100
    sd.ley.low <- mean.ley-sd.ley
    sd.ley.up <- mean.ley+sd.ley
    
    dt.T <- read.table(paste0(path.main,'thresholds_temperature_smallsample_v6.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.T[,1])
    
    seq.sample <- seq(from=100,to=10,by=-10)
    
    x <- 120
    for(i in 1:length(seq.sample))
    {
      x <- c(x,seq.sample[i])
      read_ <- dt.T[,'n.samples'] == seq.sample[i]
      
      mean.kim.loop <- mean(dt.T[read_,'kim']/dt.T[read_,'n'])*100
      sd.kim.loop <- sd(dt.T[read_,'kim']/dt.T[read_,'n'])*100
      
      mean.kim <- c(mean.kim,mean.kim.loop)
      sd.kim.low <- c(sd.kim.low,mean.kim.loop-sd.kim.loop)
      sd.kim.up <- c(sd.kim.up,mean.kim.loop+sd.kim.loop)
      
      mean.hub.loop <- mean(dt.T[read_,'hub']/dt.T[read_,'n'])*100
      sd.hub.loop <- sd(dt.T[read_,'hub']/dt.T[read_,'n'])*100
      
      mean.hub <- c(mean.hub,mean.hub.loop)
      sd.hub.low <- c(sd.hub.low,mean.hub.loop-sd.hub.loop)
      sd.hub.up <- c(sd.hub.up,mean.hub.loop+sd.hub.loop)
      
      
      mean.sch.loop <- mean(dt.T[read_,'sch']/dt.T[read_,'n'])*100
      sd.sch.loop <- sd(dt.T[read_,'sch']/dt.T[read_,'n'])*100
      
      mean.sch <- c(mean.sch,mean.sch.loop)
      sd.sch.low <- c(sd.sch.low,mean.sch.loop-sd.sch.loop)
      sd.sch.up <- c(sd.sch.up,mean.sch.loop+sd.sch.loop)
      
      
      mean.ley.loop <- mean(dt.T[read_,'ley']/dt.T[read_,'n'])*100
      sd.ley.loop <- sd(dt.T[read_,'ley']/dt.T[read_,'n'])*100
      
      mean.ley <- c(mean.ley,mean.ley.loop)
      sd.ley.low <- c(sd.ley.low,mean.ley.loop-sd.ley.loop)
      sd.ley.up <- c(sd.ley.up,mean.ley.loop+sd.ley.loop)
      
      
    }
    
    data.out <- data.frame(sample.size=x,sd.kim.low=sd.kim.low,mean.kim=mean.kim,sd.kim.up=sd.kim.up,sd.sch.low=sd.sch.low,mean.sch=mean.sch,sd.sch.up=sd.sch.up,sd.ley.low=sd.ley.low,mean.ley=mean.ley,sd.ley.up=sd.ley.up,sd.hub.low=sd.hub.low,mean.hub=mean.hub,sd.hub.up=sd.hub.up,stringsAsFactors = F)
    data.out <- data.out[order(data.out[,'sample.size']),]
    
    
    plot(data.out[,'sample.size'],data.out[,'mean.kim'],ylim=c(0,10),col=colkim,type='l',ylab='temp')
    lines(data.out[,'sample.size'],data.out[,'sd.kim.low'],col=colkim,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.kim.up'],col=colkim,lty=2)
    
    lines(data.out[,'sample.size'],data.out[,'mean.hub'],col=colhub)
    lines(data.out[,'sample.size'],data.out[,'sd.hub.low'],col=colhub,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.hub.up'],col=colhub,lty=2)
    
    lines(data.out[,'sample.size'],data.out[,'mean.sch'],col=colsch)
    lines(data.out[,'sample.size'],data.out[,'sd.sch.low'],col=colsch,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.sch.up'],col=colsch,lty=2)
    
    lines(data.out[,'sample.size'],data.out[,'mean.ley'],col=colley)
    lines(data.out[,'sample.size'],data.out[,'sd.ley.low'],col=colley,lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.ley.up'],col=colley,lty=2)
  }
}


if(1) # plots for methods logbox
{
  if(1) # precipitation
  {
    dt.P <- read.table(paste0(path.main,'thresholds_precipitation_v4.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.P[,1])
    
    # subsampling methods : sqrt(n) blocks subsampled sqrt(n) times without replacement
    N.blocks <- round(sqrt(n),digits=0)
    N.loop <- N.blocks
    list.index <- list()
    all.index <- 1:n
    for(i in 1:(N.loop-1))
    {
      subindex <- sample(all.index,N.blocks,replace = FALSE)
      all.index <- all.index[!(all.index %in% subindex)]
      list.index[[i]] <- subindex
    }
    list.index[[N.loop]] <- all.index
    
    all.mean.P <- c()
    true.percentage.threshold <- c()
    for(i in 1:N.loop)
    {
      dt.P.loop <- dt.P[list.index[[i]],]
      all.mean.P <- c(all.mean.P,(sum(dt.P.loop[,'logbox'])/sum(dt.P.loop[,'n']))*100)
      true.percentage.threshold <- c(true.percentage.threshold,median(0.001*sqrt(dt.P.loop[,'n'])*100/dt.P.loop[,'n']))
    }
    
    
    mean.logbox.subsample <- mean(all.mean.P)
    mean.logbox.real <- (sum(dt.P[,'logbox'])/sum(dt.P[,'n']))*100
    sd.logbox <- sd(all.mean.P)
    sd.logbox.low.P <- mean.logbox.real-sd.logbox
    sd.logbox.up.P <- mean.logbox.real+sd.logbox
    mean.logbox.P <- mean.logbox.real
    
    median.threshold.P <- median(true.percentage.threshold)

    dt.P <- read.table(paste0(path.main,'thresholds_precipitation_smallsample_v6.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.P[dt.P[,'n.samples'] == 10,1])
    
    # subsampling methods : sqrt(n) blocks subsampled sqrt(n) times without replacement
    N.blocks <- round(sqrt(n),digits=0)
    N.loop <- N.blocks
    list.index <- list()
    all.index <- 1:n
    for(i in 1:(N.loop-1))
    {
      subindex <- sample(all.index,N.blocks,replace = FALSE)
      all.index <- all.index[!(all.index %in% subindex)]
      list.index[[i]] <- subindex
    }
    list.index[[N.loop]] <- all.index
    
    
    
    seq.sample <- seq(from=100,to=10,by=-10)
    
    x <- c()
    sd.logbox.low <- c()
    sd.logbox.up <- c()
    mean.logbox <- c()
    for(i in 1:length(seq.sample))
    {
      x <- c(x,seq.sample[i])
      read_ <- dt.P[,'n.samples'] == seq.sample[i]
      
      dt.temp <- dt.P[read_,]
      
      all.mean <- c()
      for(j in 1:N.loop)
      {
        dt.P.loop <- dt.temp[list.index[[j]],]
        all.mean <- c(all.mean,(sum(dt.P.loop[,'logbox'])/sum(dt.P.loop[,'n']))*100)
      }
      
      sd.logbox <- sd(all.mean)
      mean.logbox.real <- (sum(dt.temp[,'logbox'])/sum(dt.temp[,'n']))*100
      
      sd.logbox.low <- c(sd.logbox.low,mean.logbox.real-sd.logbox)
      sd.logbox.up <- c(sd.logbox.up,mean.logbox.real+sd.logbox)
      mean.logbox <- c(mean.logbox,mean.logbox.real)
    }
    
    data.out <- data.frame(sample.size=x,sd.logbox.low=sd.logbox.low,mean.logbox=mean.logbox,sd.logbox.up=sd.logbox.up,stringsAsFactors = F)
    data.out <- data.out[order(data.out[,'sample.size']),]
    
    plot(data.out[,'sample.size'],data.out[,'mean.logbox'],col='blue',type='l',ylim=c(0,0.3))
    lines(data.out[,'sample.size'],data.out[,'sd.logbox.low'],col='blue',lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.logbox.up'],col='blue',lty=2)
    
    mean.ok <- 0.001*sqrt(data.out[,'sample.size'])*100/data.out[,'sample.size']
    
    lines(data.out[,'sample.size'],mean.ok,col='green')
    
  }
  
  if(1) # temperature
  {
    dt.T <- read.table(paste0(path.main,'thresholds_temperature_v4.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.T[,1])
    
    # subsampling methods : sqrt(n) blocks subsampled sqrt(n) times without replacement
    N.blocks <- round(sqrt(n),digits=0)
    N.loop <- N.blocks
    list.index <- list()
    all.index <- 1:n
    for(i in 1:(N.loop-1))
    {
      subindex <- sample(all.index,N.blocks,replace = FALSE)
      all.index <- all.index[!(all.index %in% subindex)]
      list.index[[i]] <- subindex
    }
    list.index[[N.loop]] <- all.index
    
    all.mean.T <- c()
    true.percentage.threshold <- c()
    for(i in 1:N.loop)
    {
      dt.T.loop <- dt.T[list.index[[i]],]
      all.mean.T <- c(all.mean.T,(sum(dt.T.loop[,'logbox'])/sum(dt.T.loop[,'n']))*100)
      true.percentage.threshold <- c(true.percentage.threshold,median(0.001*sqrt(dt.T.loop[,'n'])*100/dt.T.loop[,'n']))
    }
    
    mean.logbox.subsample <- mean(all.mean.T)
    mean.logbox.real <- (sum(dt.T[,'logbox'])/sum(dt.T[,'n']))*100
    sd.logbox <- sd(all.mean.T)
    sd.logbox.low.T <- mean.logbox.real-sd.logbox
    sd.logbox.up.T <- mean.logbox.real+sd.logbox
    mean.logbox.T <- mean.logbox.real
    
    median.threshold.T <- median(true.percentage.threshold)

    
    dt.T <- read.table(paste0(path.main,'thresholds_temperature_smallsample_v6.csv'),sep=',',stringsAsFactors = F,header=T)
    n <- length(dt.T[dt.T[,'n.samples'] == 10,1])
    
    # subsampling methods : sqrt(n) blocks subsampled sqrt(n) times without replacement
    N.blocks <- round(sqrt(n),digits=0)
    N.loop <- N.blocks
    list.index <- list()
    all.index <- 1:n
    for(i in 1:(N.loop-1))
    {
      subindex <- sample(all.index,N.blocks,replace = FALSE)
      all.index <- all.index[!(all.index %in% subindex)]
      list.index[[i]] <- subindex
    }
    list.index[[N.loop]] <- all.index
    
    
    
    seq.sample <- seq(from=100,to=10,by=-10)
    
    x <- c()
    sd.logbox.low <- c()
    sd.logbox.up <- c()
    mean.logbox <- c()
    for(i in 1:length(seq.sample))
    {
      x <- c(x,seq.sample[i])
      read_ <- dt.T[,'n.samples'] == seq.sample[i]
      
      dt.temp <- dt.T[read_,]
      
      all.mean <- c()
      for(j in 1:N.loop)
      {
        dt.T.loop <- dt.temp[list.index[[j]],]
        all.mean <- c(all.mean,(sum(dt.T.loop[,'logbox'])/sum(dt.T.loop[,'n']))*100)
      }
      
      sd.logbox <- sd(all.mean)
      mean.logbox.real <- (sum(dt.temp[,'logbox'])/sum(dt.temp[,'n']))*100
      
      sd.logbox.low <- c(sd.logbox.low,mean.logbox.real-sd.logbox)
      sd.logbox.up <- c(sd.logbox.up,mean.logbox.real+sd.logbox)
      mean.logbox <- c(mean.logbox,mean.logbox.real)
    }
    
    data.out <- data.frame(sample.size=x,sd.logbox.low=sd.logbox.low,mean.logbox=mean.logbox,sd.logbox.up=sd.logbox.up,stringsAsFactors = F)
    data.out <- data.out[order(data.out[,'sample.size']),]
    
    lines(data.out[,'sample.size'],data.out[,'mean.logbox'],col='red')
    lines(data.out[,'sample.size'],data.out[,'sd.logbox.low'],col='red',lty=2)
    lines(data.out[,'sample.size'],data.out[,'sd.logbox.up'],col='red',lty=2)
    

    
  }
  
  
  q.0.025.P <- quantile(all.mean.P,0.025)
  q.0.25.P <- quantile(all.mean.P,0.25)
  q.0.50.P <- quantile(all.mean.P,0.5)
  q.0.75.P <- quantile(all.mean.P,0.75)
  q.0.975.P <- quantile(all.mean.P,0.975)
  
  y.P <- c(q.0.025.P,q.0.25.P,q.0.50.P,q.0.75.P,q.0.975.P)
  
  
  q.0.025.T <- quantile(all.mean.T,0.025)
  q.0.25.T <- quantile(all.mean.T,0.25)
  q.0.50.T <- quantile(all.mean.T,0.5)
  q.0.75.T <- quantile(all.mean.T,0.75)
  q.0.975.T <- quantile(all.mean.T,0.975)
  
  y.T <- c(q.0.025.T,q.0.25.T,q.0.50.T,q.0.75.T,q.0.975.T)
  
  
  plot(rep(1,5),y.P,xlim=c(1,2),ylim=c(0,0.006),col="blue")
  points(rep(2,5),y.T,col='red')
  abline(h=median.threshold.P,col='blue')
  abline(h=median.threshold.T,col='red')

}








