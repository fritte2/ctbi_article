# This scripts calculates for each distribution part of the pearson family the alpha values, skewness and kurtosis.
# Adequate shapes values were approximated beforehand using a MonteCarlo simulation.

# working folder
setwd("~/Science/hess-2021-609/PartI.outliers/step1_alpha_calculation")

library(gsl) # required for pearsonDS
library(PearsonDS) # for the pearson distribution
library(stats) # for the gamma and beta distribution
library(extraDistr) # for the inverse gamma and betaprime distribution
library(mrfDepth) # for the med couple

# clarify which package to use
rinvgamma <- extraDistr::rinvgamma
qinvgamma <- extraDistr::qinvgamma

p.lim.3sigma = pnorm(3)
p.lim.4sigma = pnorm(4)
p.lim.5sigma = pnorm(5)

N.size <- 10^6

data <- read.table('data.grid_v9.csv',header=T,sep=',',stringsAsFactors = F)
n <- length(data[,1])
data <- data.frame(data,MC=rep(NA,n),alpha.minus.3sigma=rep(NA,n),alpha.plus.3sigma=rep(NA,n),alpha.minus.4sigma=rep(NA,n),alpha.plus.4sigma=rep(NA,n),alpha.minus.5sigma=rep(NA,n),alpha.plus.5sigma=rep(NA,n),stringsAsFactors = F)

# Functions to calculate the kurtosis, skewness, alpha.minus and alpha.plus of the distributions
log.function <- 1
if(log.function)
{
  # normal
  get.kurtosis.skewness.alpha.normal <- function(plim)
  {
    q0.75 <- qnorm(0.75)
    q0.25 <- qnorm(0.25)
    qplus <- qnorm(plim)
    qminus <- qnorm(1-plim)
    
    kur <- 3
    
    skew <- 0
    
    alpha.minus <- (q0.25-qminus)/(q0.75-q0.25)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    return(c(kur,skew,alpha.minus,alpha.plus))
  }
  
  # student
  get.kurtosis.skewness.alpha.student <- function(plim,nu)
  {
    q0.75 <- qt(0.75,df=nu)
    q0.25 <- qt(0.25,df=nu)
    qplus <- qt(plim,df=nu)
    qminus <- qt(1-plim,df=nu)
    
    kur <- 3+6/(nu-4)
    
    skew <- 0
    
    alpha.minus <- (q0.25-qminus)/(q0.75-q0.25)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    return(c(kur,skew,alpha.minus,alpha.plus))
  }
  
  # pearson IV
  get.kurtosis.skewness.alpha.pearsonIV <- function(plim,k1,k2)
  {
    q0.75 <- qpearsonIV(0.75, m=k1, nu=k2,location=0,scale=1)
    q0.25 <- qpearsonIV(0.25, m=k1, nu=k2,location=0,scale=1)
    qplus <- qpearsonIV(plim, m=k1, nu=k2,location=0,scale=1)
    qminus <- qpearsonIV(1-plim, m=k1, nu=k2,location=0,scale=1)
    
    m <- k1
    nu <- k2
    r0 <- 2*(m-1)
    
    kur <- 3*(r0-1)*((r0+6)*(r0^2+nu^2)-8*(r0^2))/((r0-2)*(r0-3)*(r0^2+nu^2))
    
    skew <- -4*nu*sqrt((r0-1)/(r0^2+nu^2))/(r0-2)
    
    alpha.minus <- (q0.25-qminus)/(q0.75-q0.25)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    return(c(kur,skew,alpha.minus,alpha.plus))
  }
  
  # gamma
  get.kurtosis.skewness.alpha.gamma <- function(plim,k)
  {
    q0.75 <- qgamma(0.75,shape=k)
    q0.25 <- qgamma(0.25,shape=k)
    qplus <- qgamma(plim,shape=k)
    qminus <- qgamma(1-plim,shape=k)
    
    kur <- 3+6/k
    
    skew <- 2/sqrt(k)
    
    alpha.minus <- (q0.25-qminus)/(q0.75-q0.25)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    return(c(kur,skew,alpha.minus,alpha.plus))
  }
  
  # invgamma
  get.kurtosis.skewness.alpha.invgamma <- function(plim,k)
  {
    q0.75 <- qinvgamma(0.75, alpha=k)
    q0.25 <- qinvgamma(0.25, alpha=k)
    qplus <- qinvgamma(plim, alpha=k)
    qminus <- qinvgamma(1-plim, alpha=k)
    
    kur <- 3+6*(5*k-11)/((k-3)*(k-4))
    
    skew <- 4*sqrt(k-2)/(k-3)
    
    alpha.minus <- (q0.25-qminus)/(q0.75-q0.25)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    return(c(kur,skew,alpha.minus,alpha.plus))
  }
  
  # beta
  get.kurtosis.skewness.alpha.beta <- function(plim,k1,k2)
  {
    q0.75 <- qbeta(0.75, shape1=k1, shape2=k2)
    q0.25 <- qbeta(0.25, shape1=k1, shape2=k2)
    qplus <- qbeta(plim, shape1=k1, shape2=k2)
    qminus <- qbeta(1-plim, shape1=k1, shape2=k2)
    
    kur <- 3+6*(((k1-k2)^2)*(k1+k2+1)-k1*k2*(k1+k2+2))/(k1*k2*(k1+k2+2)*(k1+k2+3))
    
    skew <- 2*(k2-k1)*sqrt(1+k1+k2)/((k1+k2+2)*sqrt(k1*k2))
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    alpha.minus <- (q0.25-qminus)/(q0.75-q0.25)
    
    return(c(kur,skew,alpha.minus,alpha.plus))
  }
  
  # beta prime
  get.kurtosis.skewness.alpha.betaprime <- function(plim,k1,k2)
  {
    q0.75 <- qbetapr(0.75, shape1=k1, shape2=k2)
    q0.25 <- qbetapr(0.25, shape1=k1, shape2=k2)
    qplus <- qbetapr(plim, shape1=k1, shape2=k2)
    qminus <- qbetapr(1-plim, shape1=k1, shape2=k2)
    
    kur <- 3+6*(k1*(k1+k2-1)*(5*k2-11)+((k2-1)^2)*(k2-2))/(k1*(k1+k2-1)*(k2-3)*(k2-4))
    
    skew <- 2*(2*k1+k2-1)*sqrt((k2-2)/(k1*(k1+k2-1)))/(k2-3)
    
    alpha.minus <- (q0.25-qminus)/(q0.75-q0.25)
    
    alpha.plus <- (qplus-q0.75)/(q0.75-q0.25)
    
    return(c(kur,skew,alpha.minus,alpha.plus))
  }
}


for(i in 1:n)
{
  print(paste0(round(i*100/n,digits=1),' %'))
  
  shape1 <- data[i,'shape1']
  shape2 <- data[i,'shape2']
  char.test <- data[i,'distribution']
  
  if(char.test=='normal')
  {
    d.test <- rnorm(N.size)
    ksa.3sigma <- get.kurtosis.skewness.alpha.normal(p.lim.3sigma)
    ksa.4sigma <- get.kurtosis.skewness.alpha.normal(p.lim.4sigma)
    ksa.5sigma <- get.kurtosis.skewness.alpha.normal(p.lim.5sigma)
  }
  if(char.test=='student')
  {
    d.test <- rt(N.size, df=shape1)
    ksa.3sigma <- get.kurtosis.skewness.alpha.student(p.lim.3sigma,shape1)
    ksa.4sigma <- get.kurtosis.skewness.alpha.student(p.lim.4sigma,shape1)
    ksa.5sigma <- get.kurtosis.skewness.alpha.student(p.lim.5sigma,shape1)
  }
  if(char.test=='gamma')
  {
    d.test <- rgamma(N.size,shape=shape1)
    ksa.3sigma <- get.kurtosis.skewness.alpha.gamma(p.lim.3sigma,shape1)
    ksa.4sigma <- get.kurtosis.skewness.alpha.gamma(p.lim.4sigma,shape1)
    ksa.5sigma <- get.kurtosis.skewness.alpha.gamma(p.lim.5sigma,shape1)
  }
  if(char.test=='inverse gamma')
  {
    print(shape2)
    d.test <- rinvgamma(N.size,alpha=shape1)
    ksa.3sigma <- get.kurtosis.skewness.alpha.invgamma(p.lim.3sigma,shape1)
    ksa.4sigma <- get.kurtosis.skewness.alpha.invgamma(p.lim.4sigma,shape1)
    ksa.5sigma <- get.kurtosis.skewness.alpha.invgamma(p.lim.5sigma,shape1)
  }
  if(char.test=='pearson IV')
  {
    d.test <- rpearsonIV(N.size,m=shape1,nu=shape2,location=0,scale=1)
    ksa.3sigma <- get.kurtosis.skewness.alpha.pearsonIV(p.lim.3sigma,shape1,shape2)
    ksa.4sigma <- get.kurtosis.skewness.alpha.pearsonIV(p.lim.4sigma,shape1,shape2)
    ksa.5sigma <- get.kurtosis.skewness.alpha.pearsonIV(p.lim.5sigma,shape1,shape2)
  }
  if(char.test=='beta')
  {
    d.test <- rbeta(N.size,shape1=shape1,shape2=shape2)
    ksa.3sigma <- get.kurtosis.skewness.alpha.beta(p.lim.3sigma,shape1,shape2)
    ksa.4sigma <- get.kurtosis.skewness.alpha.beta(p.lim.4sigma,shape1,shape2)
    ksa.5sigma <- get.kurtosis.skewness.alpha.beta(p.lim.5sigma,shape1,shape2)
  }
  if(char.test=='betaprime')
  {
    d.test <- rbetapr(N.size,shape1=shape1,shape2=shape2)
    ksa.3sigma <- get.kurtosis.skewness.alpha.betaprime(p.lim.3sigma,shape1,shape2)
    ksa.4sigma <- get.kurtosis.skewness.alpha.betaprime(p.lim.4sigma,shape1,shape2)
    ksa.5sigma <- get.kurtosis.skewness.alpha.betaprime(p.lim.5sigma,shape1,shape2)
  }
  
  
  data[i,'alpha.minus.3sigma'] <- round(ksa.3sigma[3],digits=4)
  data[i,'alpha.plus.3sigma'] <- round(ksa.3sigma[4],digits=4)
  data[i,'alpha.minus.4sigma'] <- round(ksa.4sigma[3],digits=4)
  data[i,'alpha.plus.4sigma'] <- round(ksa.4sigma[4],digits=4)
  data[i,'alpha.minus.5sigma'] <- round(ksa.5sigma[3],digits=4)
  data[i,'alpha.plus.5sigma'] <- round(ksa.5sigma[4],digits=4)

  MC <- as.numeric(medcouple(d.test))
  data[i,'MC'] <- round(MC,digits=4)
  
  # Analytic kurtosis & skewness
  kurtosis <- ksa.3sigma[1]
  skewness <- ksa.3sigma[2]
  data[i,'skewness'] <- round(skewness,digits=4)
  data[i,'kurtosis'] <- round(kurtosis,digits=4)
}


write.table(data,'data.alpha_vNONE.csv',sep=',',row.names = F)