# This scripts compares the 5 different models to flag outliers
# for a different sample size and different distributions.

# working folder
setwd("~/Science/hess-2021-609/PartI.outliers/step2_comparison")

source('f.Kimber_v1.R')
source('f.Hubert_v1.R')
source('f.LogBox_v1.R')
source('f.Leys_v1.R')
source('f.Schwertman_v1.R')

library(gsl) # for the pearson distribution
library(PearsonDS) # for the pearson distribution
library(stats) # for the student, gamma and beta distribution
library(extraDistr) # for the inverse gamma and betaprime distribution
library(mrfDepth) # for the med couple

data <- read.table('data.alpha_v2.csv',sep=',',header=T,stringsAsFactors = F)

sample.size <- 2^(4:14)
N.realization <- 10^6
Z.schwertman <- 3
k.coeff.outliers <- 0.60

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

# pick 100 random distributions from the student, gamma, inverse-gamma, beta, betaprime and pearson
# This is done to harmonize weights
set.seed(1)
N.100 <- 100
data.student <- data['student'==data[,1],]
data.gamma <- data['gamma'==data[,1],]
data.inversegamma <- data['inverse gamma'==data[,1],]
data.beta <- data['beta'==data[,1],]
data.betaprime <- data['betaprime'==data[,1],]
data.pearson <- data['pearson IV'==data[,1],]

data.student <- data.student[sample(1:length(data.student[,1]),N.100),]
data.gamma <- data.gamma[sample(1:length(data.gamma[,1]),N.100),]
data.inversegamma <- data.inversegamma[sample(1:length(data.inversegamma[,1]),N.100),]
data.beta <- data.beta[sample(1:length(data.beta[,1]),N.100),]
data.betaprime <- data.betaprime[sample(1:length(data.betaprime[,1]),N.100),]
data.pearson <- data.pearson[sample(1:length(data.pearson[,1]),N.100),]

data <- rbind(data.student,data.gamma,data.inversegamma,data.beta,data.betaprime,data.pearson)
n.dist <- length(data[,1])

write.table(data,'data.alpha.subset.600_v1.csv',sep=',',row.names = F)

# Store the 0.025,0.25,0.5,0.75 and 0.975 quantiles
data.0.025 <- data.frame(models=c('Kimber','Hubert','Schwertman','Leys','LogBox'),Nsize.16=rep(NA,5),Nsize.32=rep(NA,5),Nsize.64=rep(NA,5),Nsize.128=rep(NA,5),Nsize.256=rep(NA,5),Nsize.512=rep(NA,5),Nsize.1024=rep(NA,5),Nsize.2048=rep(NA,5),Nsize.4096=rep(NA,5),Nsize.8192=rep(NA,5),Nsize.16384=rep(NA,5),stringsAsFactors = F)
data.0.25 <- data.frame(models=c('Kimber','Hubert','Schwertman','Leys','LogBox'),Nsize.16=rep(NA,5),Nsize.32=rep(NA,5),Nsize.64=rep(NA,5),Nsize.128=rep(NA,5),Nsize.256=rep(NA,5),Nsize.512=rep(NA,5),Nsize.1024=rep(NA,5),Nsize.2048=rep(NA,5),Nsize.4096=rep(NA,5),Nsize.8192=rep(NA,5),Nsize.16384=rep(NA,5),stringsAsFactors = F)
data.0.50 <- data.frame(models=c('Kimber','Hubert','Schwertman','Leys','LogBox'),Nsize.16=rep(NA,5),Nsize.32=rep(NA,5),Nsize.64=rep(NA,5),Nsize.128=rep(NA,5),Nsize.256=rep(NA,5),Nsize.512=rep(NA,5),Nsize.1024=rep(NA,5),Nsize.2048=rep(NA,5),Nsize.4096=rep(NA,5),Nsize.8192=rep(NA,5),Nsize.16384=rep(NA,5),stringsAsFactors = F)
data.0.75 <- data.frame(models=c('Kimber','Hubert','Schwertman','Leys','LogBox'),Nsize.16=rep(NA,5),Nsize.32=rep(NA,5),Nsize.64=rep(NA,5),Nsize.128=rep(NA,5),Nsize.256=rep(NA,5),Nsize.512=rep(NA,5),Nsize.1024=rep(NA,5),Nsize.2048=rep(NA,5),Nsize.4096=rep(NA,5),Nsize.8192=rep(NA,5),Nsize.16384=rep(NA,5),stringsAsFactors = F)
data.0.975 <- data.frame(models=c('Kimber','Hubert','Schwertman','Leys','LogBox'),Nsize.16=rep(NA,5),Nsize.32=rep(NA,5),Nsize.64=rep(NA,5),Nsize.128=rep(NA,5),Nsize.256=rep(NA,5),Nsize.512=rep(NA,5),Nsize.1024=rep(NA,5),Nsize.2048=rep(NA,5),Nsize.4096=rep(NA,5),Nsize.8192=rep(NA,5),Nsize.16384=rep(NA,5),stringsAsFactors = F)



for(k.size in 1:length(sample.size))
{
  n.size <- sample.size[k.size]
  
  per.Kimber <- NULL
  per.Hubert <- NULL
  per.Schwertman <- NULL
  per.Leys <- NULL
  per.LogBox <- NULL
  
for(k.dist in 1:n.dist)
{
  print(paste0(round(k.dist*100/n.dist,digits=1),' %'))
  
  shape1.test <- data[k.dist,'shape1']
  shape2.test <- data[k.dist,'shape2']
  char.test <- data[k.dist,'distribution']

  N.rep <- ceiling(N.realization/n.size)
  
  n.all <- 0

  ncaptured.Kimber <- 0
  ncaptured.Hubert <- 0
  ncaptured.Schwertman <- 0
  ncaptured.Leys <- 0
  ncaptured.LogBox <- 0
  
  for(k.rep in 1:N.rep)
  {
    
  if(char.test=='student')
  {
    d.loop <- rt(n.size, df=shape1.test)
  }
  if(char.test=='gamma')
  {
    d.loop <- rgamma(n.size,shape=shape1.test)
  }
  if(char.test=='inverse gamma')
  {
    d.loop <- rinvgamma(n.size,alpha=shape1.test)
  }
  if(char.test=='pearson IV')
  {
    d.loop <- rpearsonIV(n.size,m=shape1.test,nu=shape2.test,location=0,scale=1)
  }
  if(char.test=='beta')
  {
    d.loop <- rbeta(n.size,shape1=shape1.test,shape2=shape2.test)
  }
  if(char.test=='betaprime')
  {
    d.loop <- rbetapr(n.size,shape1=shape1.test,shape2=shape2.test)
  }

  n.all <- n.all+n.size
  MC <- as.numeric(medcouple(d.loop))
  MAD0 <- as.numeric(mad(d.loop))
  q0.25 <- as.numeric(quantile(d.loop,0.25))
  q0.50 <- as.numeric(quantile(d.loop,0.50))
  q0.75 <- as.numeric(quantile(d.loop,0.75))
  
  ncaptured.Kimber.loop <- f.Kimber(d.loop,q0.25,q0.50,q0.75)
  ncaptured.Hubert.loop <- f.Hubert(d.loop,q0.25,q0.75,MC)
  ncaptured.Schwertman.loop <- f.Schwertman(d.loop,q0.25,q0.50,q0.75,Z.schwertman,k.n.Schwertman[k.size])
  ncaptured.Leys.loop <- f.Leys(d.loop,q0.50,MAD0)
  ncaptured.LogBox.loop <- f.LogBox(d.loop,q0.25,q0.75,k.coeff.outliers)
  
  ncaptured.Kimber <- ncaptured.Kimber+ncaptured.Kimber.loop
  ncaptured.Hubert <- ncaptured.Hubert+ncaptured.Hubert.loop
  ncaptured.Schwertman <- ncaptured.Schwertman+ncaptured.Schwertman.loop
  ncaptured.Leys <- ncaptured.Leys+ncaptured.Leys.loop
  ncaptured.LogBox <- ncaptured.LogBox+ncaptured.LogBox.loop
  }
  
  per.Kimber <- c(per.Kimber,ncaptured.Kimber/n.all)
  per.Hubert <- c(per.Hubert,ncaptured.Hubert/n.all)
  per.Schwertman <- c(per.Schwertman,ncaptured.Schwertman/n.all)
  per.Leys <- c(per.Leys,ncaptured.Leys/n.all)
  per.LogBox <- c(per.LogBox,ncaptured.LogBox/n.all)
}
  
  data.0.025[1,k.size+1] <- quantile(per.Kimber,0.025)
  data.0.025[2,k.size+1] <- quantile(per.Hubert,0.025)
  data.0.025[3,k.size+1] <- quantile(per.Schwertman,0.025)
  data.0.025[4,k.size+1] <- quantile(per.Leys,0.025)
  data.0.025[5,k.size+1] <- quantile(per.LogBox,0.025)
  
  data.0.25[1,k.size+1] <- quantile(per.Kimber,0.25)
  data.0.25[2,k.size+1] <- quantile(per.Hubert,0.25)
  data.0.25[3,k.size+1] <- quantile(per.Schwertman,0.25)
  data.0.25[4,k.size+1] <- quantile(per.Leys,0.25)
  data.0.25[5,k.size+1] <- quantile(per.LogBox,0.25)
  
  data.0.50[1,k.size+1] <- quantile(per.Kimber,0.5)
  data.0.50[2,k.size+1] <- quantile(per.Hubert,0.5)
  data.0.50[3,k.size+1] <- quantile(per.Schwertman,0.5)
  data.0.50[4,k.size+1] <- quantile(per.Leys,0.5)
  data.0.50[5,k.size+1] <- quantile(per.LogBox,0.5)
  
  data.0.75[1,k.size+1] <- quantile(per.Kimber,0.75)
  data.0.75[2,k.size+1] <- quantile(per.Hubert,0.75)
  data.0.75[3,k.size+1] <- quantile(per.Schwertman,0.75)
  data.0.75[4,k.size+1] <- quantile(per.Leys,0.75)
  data.0.75[5,k.size+1] <- quantile(per.LogBox,0.75)
  
  data.0.975[1,k.size+1] <- quantile(per.Kimber,0.975)
  data.0.975[2,k.size+1] <- quantile(per.Hubert,0.975)
  data.0.975[3,k.size+1] <- quantile(per.Schwertman,0.975)
  data.0.975[4,k.size+1] <- quantile(per.Leys,0.975)
  data.0.975[5,k.size+1] <- quantile(per.LogBox,0.975)
}

write.table(data.0.025,'data.0.025_v3.csv',row.names = F,sep=',')
write.table(data.0.25,'data.0.25_v3.csv',row.names = F,sep=',')
write.table(data.0.50,'data.0.50_v3.csv',row.names = F,sep=',')
write.table(data.0.75,'data.0.75_v3.csv',row.names = F,sep=',')
write.table(data.0.975,'data.0.975_v3.csv',row.names = F,sep=',')





