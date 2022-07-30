# working folder
setwd("~/Science/ctbi12/PartI.outliers/step6_plots/FIG1/")

spatial.agg <- function(df,xy.step,xlim,ylim)
{
  x.start <- xlim[1]+xy.step/2
  y.start <- ylim[1]+xy.step/2
  
  seq.x <- seq(from=x.start,to=xlim[2],by=xy.step)
  seq.y <- seq(from=y.start,to=ylim[2],by=xy.step)
  N <- length(seq.x)*length(seq.y)
  df.agg <- data.frame(x=rep(NA,N),y=rep(NA,N),A=rep(NA,N),B=rep(NA,N),stringsAsFactors = F)
  k.loop <- 0
  
  for(i in 1:length(seq.x))
  {
    for(j in 1:length(seq.y))
    {
      read_ <- ((seq.x[i]-xy.step/2) <= df[,'x'] & df[,'x'] < (seq.x[i]+xy.step/2))  & ((seq.y[j]-xy.step/2) <= df[,'y'] & df[,'y'] < (seq.y[j]+xy.step/2))
      k.loop <- k.loop+1
      df.agg[k.loop,'x'] <- seq.x[i]
      df.agg[k.loop,'y'] <- seq.y[j]
      if(sum(read_) > 0)
      {
        df.agg[k.loop,'A'] <- mean(df[read_,'A'])
        df.agg[k.loop,'B'] <- mean(df[read_,'B'])
      }else
      {
        df.agg[k.loop,'A'] <- NA
        df.agg[k.loop,'B'] <- NA
      }
    }
  }
  
  df.agg <- df.agg[!is.na(df.agg[,'A']),]
  
  return(df.agg)
}

library(viridis)
library(ggplot2)

data <- read.table('data.alpha_AB_v10.csv',sep=',',header=T,stringsAsFactors = F)
pearson.family <- c('normal','gamma','inverse gamma','betaprime','student','pearson IV')
data.pearson <- data[data[,'distribution'] %in% pearson.family,]

# pearson family
if(0)
{
  spatial.agg <- function(df,xy.step,xlim,ylim)
  {
    x.start <- xlim[1]+xy.step/2
    y.start <- ylim[1]+xy.step/2
    
    seq.x <- seq(from=x.start,to=xlim[2],by=xy.step)
    seq.y <- seq(from=y.start,to=ylim[2],by=xy.step)
    N <- length(seq.x)*length(seq.y)
    df.agg <- data.frame(x=rep(NA,N),y=rep(NA,N),A=rep(NA,N),B=rep(NA,N),stringsAsFactors = F)
    k.loop <- 0
    
    for(i in 1:length(seq.x))
    {
      for(j in 1:length(seq.y))
      {
        read_ <- ((seq.x[i]-xy.step/2) <= df[,'x'] & df[,'x'] < (seq.x[i]+xy.step/2))  & ((seq.y[j]-xy.step/2) <= df[,'y'] & df[,'y'] < (seq.y[j]+xy.step/2))
        k.loop <- k.loop+1
        df.agg[k.loop,'x'] <- seq.x[i]
        df.agg[k.loop,'y'] <- seq.y[j]
        if(sum(read_) > 0)
        {
          df.agg[k.loop,'A'] <- mean(df[read_,'A'])
          df.agg[k.loop,'B'] <- mean(df[read_,'B'])
        }else
        {
          df.agg[k.loop,'A'] <- NA
          df.agg[k.loop,'B'] <- NA
        }
      }
    }
    
    df.agg <- df.agg[!is.na(df.agg[,'A']),]
    
    return(df.agg)
  }
  
  df <- data.frame(x=data.pearson[,'kurtosis']-3,y=data.pearson[,'skewness']^2,A=data.pearson[,'A'],B=data.pearson[,'B'],stringsAsFactors = F)
  df <- spatial.agg(df,0.1,xlim=c(0,6),ylim=c(0,4))
  
  
  plot.A.pearson <- ggplot(df, aes(x,y)) + 
  geom_point(aes(colour = A)) + 
  scale_color_viridis(option = "D") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  

  
  plot.B.pearson <- ggplot(df,aes(x,y)) + 
  geom_point(aes(colour = B)) + 
  scale_color_viridis(option = "D") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


  
  print(plot.A.pearson)
  print(plot.B.pearson)
}

# all distributions
if(1)
{
  # there are too many points (~6000). I need to reduce the number to less than 1000.
  # I keep all the weibull points for m.star > 0.45
  # I randomly select points for each distribution below
  
  m.threshold <- 0.45
  
  data.upset <- data[data[,'m.star'] > m.threshold,]
  
  data.subset <- data[data[,'m.star'] < m.threshold,]
  
  data.agg <- data.subset[1,]
  data.agg <- data.agg[-1,]
  
  all.dist <- unique(data.subset[,'distribution'])
  
  for(i in 1:length(all.dist))
  {
    ds <- data.subset[all.dist[i] == data.subset[,'distribution'],]
    n.s <- length(ds[,1])
    
    if(n.s < 100)
    {
      data.agg <- rbind(data.agg,ds)
    }else
    {
      data.agg <- rbind(data.agg,ds[sample(1:n.s,100),])
    }
  }
  
  data.agg <- rbind(data.agg,data.upset)
  
  data.agg[data.agg[,'distribution'] != 'frechet','distribution'] <- 'good'


  x <- data.agg[,'m.star']
  a1 <- 0.2294
  a2 <- 2.9416
  a3 <- -0.0512
  a4 <- -0.0684
  y.A <- a1*exp(a2*x+a3*x^2+a4*x^3)
  
  b1 <- 1.0585
  b2 <- 15.6960
  b3 <- -17.3618
  b4 <- 28.3511
  b5 <- -11.4726
  y.B <- b1+b2*x+b3*(x^2)+b4*(x^3)+b5*(x^4)
  
  data.agg <- data.frame(data.agg,A.fit = y.A, B.fit = y.B, stringsAsFactors = F)
  data.agg <- data.agg[order(data.agg[,'m.star']),]
  
  read_ <- data.agg[,'distribution'] == 'good'
  
  plot(data.agg[read_,'m.star'],data.agg[read_,'A'],pch=16)
  lines(data.agg[read_,'m.star'],data.agg[read_,'A.fit'],col='red',lwd=2)
  points(data.agg[!read_,'m.star'],data.agg[!read_,'A'],col='skyblue',pch=16)
  
  plot(data.agg[read_,'m.star'],data.agg[read_,'B'],ylim=c(-10,25),pch=16)
  lines(data.agg[read_,'m.star'],data.agg[read_,'B.fit'],col='red',lwd=2)
  points(data.agg[!read_,'m.star'],data.agg[!read_,'B'],col='skyblue',pch=16)
  
}
