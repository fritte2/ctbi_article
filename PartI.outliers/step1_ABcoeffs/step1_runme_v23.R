# This scripts calculates for each distribution part of the Pearson and GEV family the alpha values, skewness, kurtosis, m.star, A, B and r squared of the relationship alpha = A*log(n)+B
# Adequate shapes values were approximated beforehand using a Monte-Carlo simulation (data.grid_v16.csv).

# working folder
setwd("~/Science/ctbi12/PartI.outliers/step1_ABcoeffs/")

library(gsl) # required for pearsonDS
library(PearsonDS) # for the pearson distribution
library(stats) # for the gamma distribution
library(extraDistr) # for the inverse gamma and betaprime distribution
library(VGAM) # for the weibull, frechet and gumbel distributions

# seed
set.seed(1)

# clarify which package to use
rinvgamma <- extraDistr::rinvgamma
qinvgamma <- extraDistr::qinvgamma
rfrechet <- VGAM::rfrechet
qfrechet <- VGAM::qfrechet
rgumbel <- VGAM::rgumbel
qgumbel <- VGAM::qgumbel

# alpha is calculated as alpha =  (Q(p(n))-Q(0.75))/(Q(0.75)-Q(0.25)) for each distribution (all right-skewed)
# with p(n) = 1-f(n)/(2*n) with f(n) = 0.001*sqrt(n) is the number of erroneously flagged outliers.
f.p <- function(x)
{
  p.out <- 1-0.001*sqrt(x)/(2*x)
  return(p.out)
}
p1 = f.p(100)
p2 = f.p(1000)
p3 = f.p(10000)
p4 = f.p(100000)
p5 = f.p(1000000)

N.size <- 10^6

data <- read.table('data.grid_v16.csv',header=T,sep=',',stringsAsFactors = F)
n <- length(data[,1])
data <- data.frame(data,m.star.sample=rep(NA,n),kurtosis.sample=rep(NA,n),skewness.sample=rep(NA,n),A=rep(NA,n),B=rep(NA,n),m.star=rep(NA,n),r2=rep(NA,n),alpha.100=rep(NA,n),alpha.1000=rep(NA,n),alpha.10000=rep(NA,n),alpha.100000=rep(NA,n),alpha.1000000=rep(NA,n),stringsAsFactors = F)


# Functions to calculate the kurtosis, skewness, alpha.minus and alpha and m.star of the distributions (Pearson + GEV family)
log.function <- 1
if(log.function)
{
  # moors predictor for only the longest tail (equal to half the centered moors for symmetrical distributions)
  centered.one.sided.moors <- function(res)
  {
    E1 <- quantile(res,1/8)
    E2 <- quantile(res,2/8)
    E3 <- quantile(res,3/8)
    E5 <- quantile(res,5/8)
    E6 <- quantile(res,6/8)
    E7 <- quantile(res,7/8)
    
    m1 <- (E7-E5)/(E6-E2)
    
    m2 <- (E3-E1)/(E6-E2)
    
    km <- round(max(c(m1,m2))-0.6165,digits=4)
    
    return(km)
  }
  
  # normal
  get.kurtosis.skewness.alpha.normal <- function(plim)
  {
    q0.875 <- qnorm(0.875)
    q0.75 <- qnorm(0.75)
    q0.625 <- qnorm(0.625)
    q0.25 <- qnorm(0.25)
    qplus <- qnorm(plim)
    
    kur <- 3
    
    skew <- 0
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # student
  get.kurtosis.skewness.alpha.student <- function(plim,nu)
  {
    q0.875 <- qt(0.875,df=nu)
    q0.75 <- qt(0.75,df=nu)
    q0.625 <- qt(0.625,df=nu)
    q0.25 <- qt(0.25,df=nu)
    qplus <- qt(plim,df=nu)
    
    kur <- 3+6/(nu-4)
    
    skew <- 0
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # pearson IV
  get.kurtosis.skewness.alpha.pearsonIV <- function(plim,k1,k2)
  {
    q0.875 <- qpearsonIV(0.875, m=k1, nu=k2,location=0,scale=1)
    q0.75 <- qpearsonIV(0.75, m=k1, nu=k2,location=0,scale=1)
    q0.625 <- qpearsonIV(0.625, m=k1, nu=k2,location=0,scale=1)
    q0.25 <- qpearsonIV(0.25, m=k1, nu=k2,location=0,scale=1)
    qplus <- qpearsonIV(plim, m=k1, nu=k2,location=0,scale=1)
    
    m <- k1
    nu <- k2
    r0 <- 2*(m-1)
    
    kur <- 3*(r0-1)*((r0+6)*(r0^2+nu^2)-8*(r0^2))/((r0-2)*(r0-3)*(r0^2+nu^2))
    
    skew <- -4*nu*sqrt((r0-1)/(r0^2+nu^2))/(r0-2)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # gamma
  get.kurtosis.skewness.alpha.gamma <- function(plim,k)
  {
    q0.875 <- qgamma(0.875,shape=k)
    q0.75 <- qgamma(0.75,shape=k)
    q0.625 <- qgamma(0.625,shape=k)
    q0.25 <- qgamma(0.25,shape=k)
    qplus <- qgamma(plim,shape=k)
    
    kur <- 3+6/k
    
    skew <- 2/sqrt(k)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # invgamma
  get.kurtosis.skewness.alpha.invgamma <- function(plim,k)
  {
    q0.875 <- qinvgamma(0.875, alpha=k)
    q0.75 <- qinvgamma(0.75, alpha=k)
    q0.625 <- qinvgamma(0.625, alpha=k)
    q0.25 <- qinvgamma(0.25, alpha=k)
    qplus <- qinvgamma(plim, alpha=k)
    
    kur <- 3+6*(5*k-11)/((k-3)*(k-4))
    
    skew <- 4*sqrt(k-2)/(k-3)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # beta prime
  get.kurtosis.skewness.alpha.betaprime <- function(plim,k1,k2)
  {
    q0.875 <- qbetapr(0.875, shape1=k1, shape2=k2)
    q0.75 <- qbetapr(0.75, shape1=k1, shape2=k2)
    q0.625 <- qbetapr(0.625, shape1=k1, shape2=k2)
    q0.25 <- qbetapr(0.25, shape1=k1, shape2=k2)
    qplus <- qbetapr(plim, shape1=k1, shape2=k2)
    
    kur <- 3+6*(k1*(k1+k2-1)*(5*k2-11)+((k2-1)^2)*(k2-2))/(k1*(k1+k2-1)*(k2-3)*(k2-4))
    
    skew <- 2*(2*k1+k2-1)*sqrt((k2-2)/(k1*(k1+k2-1)))/(k2-3)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # weibull
  get.kurtosis.skewness.alpha.weibull <- function(plim,k1)
  {
    q0.875 <- qweibull(0.875,shape=k1)
    q0.75 <- qweibull(0.75,shape=k1)
    q0.625 <- qweibull(0.625,shape=k1)
    q0.25 <- qweibull(0.25,shape=k1)
    qplus <- qweibull(plim,shape=k1)
    
    g1 <- gamma(1+1/k1)
    g2 <- gamma(1+2/k1)
    g3 <- gamma(1+3/k1)
    g4 <- gamma(1+4/k1)
    
    kur <- ((-6*(g1^4)+12*(g1^2)*g2-3*(g2^2)-4*g1*g3+g4)/((g2-g1^2)^2))+3
    
    skew <- (2*(g1^3)-3*g1*g2+g3)/((g2-g1^2)^(3/2))
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # gumbel
  get.kurtosis.skewness.alpha.gumbel <- function(plim)
  {
    q0.875 <- qgumbel(0.875)
    q0.75 <- qgumbel(0.75)
    q0.625 <- qgumbel(0.625)
    q0.25 <- qgumbel(0.25)
    qplus <- qgumbel(plim)
    
    kur <- 5.4
    
    skew <- 1.14
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
  
  # frechet
  get.kurtosis.skewness.alpha.frechet <- function(plim,k1)
  {
    q0.875 <- qfrechet(0.875,shape=k1)
    q0.75 <- qfrechet(0.75,shape=k1)
    q0.625 <- qfrechet(0.625,shape=k1)
    q0.25 <- qfrechet(0.25,shape=k1)
    qplus <- qfrechet(plim,shape=k1)
    
    g1 <- gamma(1-1/k1)
    g2 <- gamma(1-2/k1)
    g3 <- gamma(1-3/k1)
    g4 <- gamma(1-4/k1)
    
    kur <- ((g4-4*g3*g1+3*(g2^2))/((g2-g1^2)^2))-3
    
    skew <- (g3-3*g2*g1+2*(g1^3))/((g2-g1^2)^(3/2))
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    m.star <- (q0.875-q0.625)/(q0.75-q0.25)-0.6165
    
    return(c(kur,skew,alpha.plus,m.star))
  }
}

# produce random samples of size 10^6 just to check that the formula of kurtosis/skewness/m.star matches the observed.
for(i in 1:n)
{
  print(paste0(round(i*100/n,digits=1),' %'))
  
  shape1 <- data[i,'shape1']
  shape2 <- data[i,'shape2']
  char.test <- data[i,'distribution']
  
  if(char.test=='normal')
  {
    d.test <- rnorm(N.size)
    ksa1 <- get.kurtosis.skewness.alpha.normal(p1)
    ksa2 <- get.kurtosis.skewness.alpha.normal(p2)
    ksa3 <- get.kurtosis.skewness.alpha.normal(p3)
    ksa4 <- get.kurtosis.skewness.alpha.normal(p4)
    ksa5 <- get.kurtosis.skewness.alpha.normal(p5)
  }
  if(char.test=='student')
  {
    d.test <- rt(N.size, df=shape1)
    ksa1 <- get.kurtosis.skewness.alpha.student(p1,shape1)
    ksa2 <- get.kurtosis.skewness.alpha.student(p2,shape1)
    ksa3 <- get.kurtosis.skewness.alpha.student(p3,shape1)
    ksa4 <- get.kurtosis.skewness.alpha.student(p4,shape1)
    ksa5 <- get.kurtosis.skewness.alpha.student(p5,shape1)
  }
  if(char.test=='gamma')
  {
    d.test <- rgamma(N.size,shape=shape1)
    ksa1 <- get.kurtosis.skewness.alpha.gamma(p1,shape1)
    ksa2 <- get.kurtosis.skewness.alpha.gamma(p2,shape1)
    ksa3 <- get.kurtosis.skewness.alpha.gamma(p3,shape1)
    ksa4 <- get.kurtosis.skewness.alpha.gamma(p4,shape1)
    ksa5 <- get.kurtosis.skewness.alpha.gamma(p5,shape1)
  }
  if(char.test=='inverse gamma')
  {
    d.test <- rinvgamma(N.size,alpha=shape1)
    ksa1 <- get.kurtosis.skewness.alpha.invgamma(p1,shape1)
    ksa2 <- get.kurtosis.skewness.alpha.invgamma(p2,shape1)
    ksa3 <- get.kurtosis.skewness.alpha.invgamma(p3,shape1)
    ksa4 <- get.kurtosis.skewness.alpha.invgamma(p4,shape1)
    ksa5 <- get.kurtosis.skewness.alpha.invgamma(p5,shape1)
  }
  if(char.test=='pearson IV')
  {
    d.test <- rpearsonIV(N.size,m=shape1,nu=shape2,location=0,scale=1)
    ksa1 <- get.kurtosis.skewness.alpha.pearsonIV(p1,shape1,shape2)
    ksa2 <- get.kurtosis.skewness.alpha.pearsonIV(p2,shape1,shape2)
    ksa3 <- get.kurtosis.skewness.alpha.pearsonIV(p3,shape1,shape2)
    ksa4 <- get.kurtosis.skewness.alpha.pearsonIV(p4,shape1,shape2)
    ksa5 <- get.kurtosis.skewness.alpha.pearsonIV(p5,shape1,shape2)
  }
  if(char.test=='betaprime')
  {
    d.test <- rbetapr(N.size,shape1=shape1,shape2=shape2)
    ksa1 <- get.kurtosis.skewness.alpha.betaprime(p1,shape1,shape2)
    ksa2 <- get.kurtosis.skewness.alpha.betaprime(p2,shape1,shape2)
    ksa3 <- get.kurtosis.skewness.alpha.betaprime(p3,shape1,shape2)
    ksa4 <- get.kurtosis.skewness.alpha.betaprime(p4,shape1,shape2)
    ksa5 <- get.kurtosis.skewness.alpha.betaprime(p5,shape1,shape2)
  }
  if(char.test=='weibull')
  {
    d.test <- rweibull(N.size,shape=shape1)
    ksa1 <- get.kurtosis.skewness.alpha.weibull(p1,shape1)
    ksa2 <- get.kurtosis.skewness.alpha.weibull(p2,shape1)
    ksa3 <- get.kurtosis.skewness.alpha.weibull(p3,shape1)
    ksa4 <- get.kurtosis.skewness.alpha.weibull(p4,shape1)
    ksa5 <- get.kurtosis.skewness.alpha.weibull(p5,shape1)
  }
  if(char.test=='gumbel')
  {
    d.test <- rgumbel(N.size)
    ksa1 <- get.kurtosis.skewness.alpha.gumbel(p1)
    ksa2 <- get.kurtosis.skewness.alpha.gumbel(p2)
    ksa3 <- get.kurtosis.skewness.alpha.gumbel(p3)
    ksa4 <- get.kurtosis.skewness.alpha.gumbel(p4)
    ksa5 <- get.kurtosis.skewness.alpha.gumbel(p5)
  }
  if(char.test=='frechet')
  {
    d.test <- rfrechet(N.size,shape=shape1)
    ksa1 <- get.kurtosis.skewness.alpha.frechet(p1,shape1)
    ksa2 <- get.kurtosis.skewness.alpha.frechet(p2,shape1)
    ksa3 <- get.kurtosis.skewness.alpha.frechet(p3,shape1)
    ksa4 <- get.kurtosis.skewness.alpha.frechet(p4,shape1)
    ksa5 <- get.kurtosis.skewness.alpha.frechet(p5,shape1)
  }
  
  data[i,'alpha.100'] <- round(as.numeric(ksa1[3]),digits=4)
  data[i,'alpha.1000'] <- round(as.numeric(ksa2[3]),digits=4)
  data[i,'alpha.10000'] <- round(as.numeric(ksa3[3]),digits=4)
  data[i,'alpha.100000'] <- round(as.numeric(ksa4[3]),digits=4)
  data[i,'alpha.1000000'] <- round(as.numeric(ksa5[3]),digits=4)

  data[i,'m.star.sample'] <- round(centered.one.sided.moors(d.test),digits=4)
  data[i,'m.star'] <- ksa1[4]
  
  z.score <- (d.test-mean(d.test))/sd(d.test)
  
  skewness.sample <- mean(z.score^3)
  kurtosis.sample <- mean(z.score^4)
  
  # Analytic kurtosis & skewness
  kurtosis <- ksa1[1]
  skewness <- ksa1[2]
  data[i,'skewness'] <- round(skewness,digits=4)
  data[i,'kurtosis'] <- round(kurtosis,digits=4)
  data[i,'skewness.sample'] <- round(skewness.sample,digits=4)
  data[i,'kurtosis.sample'] <- round(kurtosis.sample,digits=4)
  
  # k and r^2 of k
  x <- log(c(100,1000,10000,100000,1000000))
  y <- as.numeric(data[i,c('alpha.100','alpha.1000','alpha.10000','alpha.100000','alpha.1000000')])
  
  k1 <- lm(y ~ x)
  
  B <- round(as.numeric(k1$coefficients[1]),digits=4)
  A <- round(as.numeric(k1$coefficients[2]),digits=4)
  

  SS.res <- sum((y-(B+A*x))^2)
  SS.tot <- sum((y-mean(y))^2)
  
  r.squared <- round(1-SS.res/SS.tot,digits=6)
  
  data[i,'A'] <- A
  data[i,'B'] <- B
  data[i,'r2'] <- r.squared
  
  if(round(i/100,digits=0) == i/100)
  {
    write.table(data,'data.alpha_AB_v10.csv',sep=',',row.names = F)
  }
  
  rm(d.test)
}

write.table(data,'data.alpha_AB_v10.csv',sep=',',row.names = F)
