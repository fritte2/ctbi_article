# this scripts calculates the optimum coefficients used to parametrize the functions g_A(m.star) and g_B(m.star) with a monte-carlo simulation.
# Fréchet is removed beforehand

setwd("~/Science/ctbi12/PartI.outliers/step3_moors/")

data <- read.table('data.alpha_AB_v10.csv',sep=',',header=T,stringsAsFactors = F)

data.nofrechet <- data[!(data[,'distribution'] %in% c('frechet')),]
data.nofrechet <- data.nofrechet[order(data.nofrechet[,'m.star']),]

# aggregate the data for m.star < 0.4 to calculate the best fit
data.subset <- data.nofrechet[data.nofrechet[,'m.star'] < 0.4,]
data.upset <- data.nofrechet[data.nofrechet[,'m.star'] > 0.4,]

m.step <- 0.025
seq.step <- seq(from=0,to=0.4,by=m.step)
n <- length(seq.step)
data.agg <- data.frame(distribution = rep('pearson',n-1),m.star=rep(NA,n-1),A=rep(NA,n-1),B=rep(NA,n-1),n.points=rep(NA,n-1),stringsAsFactors = F)
for(i in 1:(n-1))
{
  beg_ <- seq.step[i]
  end_ <- seq.step[i+1]
  read_ <- beg_ <= data.subset[,'m.star'] & data.subset[,'m.star'] < end_
  data.agg[i,'m.star'] <- median(data.subset[read_,'m.star'])
  data.agg[i,'A'] <- median(data.subset[read_,'A'])
  data.agg[i,'B'] <- median(data.subset[read_,'B'])
  data.agg[i,'n.points'] <- sum(read_)
}

data.fit <- rbind(data.agg[,c('m.star','A','B')],data.upset[,c('m.star','A','B')])

# fit for the A coefficient, monte-carlo simulation 
if(1)
{
  x <- data.fit[,'m.star']
  y <- data.fit[,'A']

  a1 <- 0.2294
  a2 <- 2.9416
  a3 <- -0.0512
  a4 <- -0.0684
  eps1 <- 0.001
  eps2 <- 0.001
  eps3 <- 0.001
  eps4 <- 0.001

  y.f <- a1*exp(a2*x+a3*x^2+a4*x^3)
  RMSE <- sqrt((1/length(y))*sum((y-y.f)^2))
  for(i in 1:1000) # can put 10^7 instead, but the coefficient are already optimized
  {
    a1.loop <- round(runif(1,min=a1-eps1,max=a1+eps1),digits=4)
    a2.loop <- round(runif(1,min=a2-eps2,max=a2+eps2),digits=4)
    a3.loop <- round(runif(1,min=a3-eps3,max=a3+eps3),digits=4)
    a4.loop <- round(runif(1,min=a4-eps4,max=a4+eps4),digits=4)

    y.f <- a1.loop*exp(a2.loop*x+a3.loop*x^2+a4.loop*x^3)
    RMSE.loop <- sqrt((1/length(y))*sum((y-y.f)^2))
    if(RMSE.loop < RMSE)
    {
      a1 <- a1.loop
      a2 <- a2.loop
      a3 <- a3.loop
      a4 <- a4.loop
      RMSE <- RMSE.loop

      print(c(a1,a2,a3,a4))
    }
  }
  plot(x,y,xlab='centered moors',ylab='A coefficient')
  y.f <- a1*exp(a2*x+a3*x^2+a4*x^3)
  lines(x,y.f,col='red',lwd=2)
  print('4 a Coefficients:')
  print(c(a1,a2,a3,a4))
  print(RMSE)

  y.f <- a1*exp(a2*x+a3*x^2+a4*x^3)
  SS.res <- sum((y-y.f)^2)
  SS.tot <- sum((y-mean(y))^2)
  r2 <- 1-SS.res/SS.tot
  print(r2)
}

# fit for the B coefficient, monte-carlo simulation 
if(1)
{
  x <- data.fit[,'m.star']
  y <- data.fit[,'B']

  b1 <- 1.0585
  b2 <- 15.6960
  b3 <- -17.3618
  b4 <- 28.3511
  b5 <- -11.4726

  eps1 <- 0.001
  eps2 <- 0.001
  eps3 <- 0.001
  eps4 <- 0.001
  eps5 <- 0.001

  y.f <- b1+b2*x+b3*(x^2)+b4*(x^3)+b5*(x^4)
  y.diff <- y-y.f
  y.diff <- y.diff[!is.na(y.diff)]
  RMSE <- sqrt((1/length(y))*sum((y.diff)^2))
  for(i in 1:1000) # can put 10^7 instead, but the coefficient are already optimized
  {
    b1.loop <- round(runif(1,min=b1-eps1,max=b1+eps1),digits=4)
    b2.loop <- round(runif(1,min=b2-eps2,max=b2+eps2),digits=4)
    b3.loop <- round(runif(1,min=b3-eps3,max=b3+eps3),digits=4)
    b4.loop <- round(runif(1,min=b4-eps4,max=b4+eps4),digits=4)
    b5.loop <- round(runif(1,min=b5-eps5,max=b5+eps5),digits=4)

    y.f <- b1.loop+b2.loop*x+b3.loop*(x^2)+b4.loop*(x^3)+b5.loop*(x^4)
    y.diff <- y-y.f
    y.diff <- y.diff[!is.na(y.diff)]
    RMSE.loop <- sqrt((1/length(y))*sum((y.diff)^2))
    if(RMSE.loop < RMSE)
    {
      b1 <- b1.loop
      b2 <- b2.loop
      b3 <- b3.loop
      b4 <- b4.loop
      b5 <- b5.loop
      RMSE <- RMSE.loop

      print(c(b1,b2,b3,b4,b5))
    }
  }
  y.f <- b1+b2*x+b3*(x^2)+b4*(x^3)+b5*(x^4)
  plot(x,y,xlab='centered moors',ylab='B coefficient')
  lines(x,y.f,col='red',lwd=2)
  print('5 b Coefficients:')
  print(c(b1,b2,b3,b4,b5))
  print(RMSE)

  SS.res <- sum((y-y.f)^2)
  SS.tot <- sum((y-mean(y))^2)
  r2 <- 1-SS.res/SS.tot
  print(r2)
}
