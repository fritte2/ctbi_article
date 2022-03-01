setwd("~/Science/hess-2021-609/PartI.outliers/step3_plots")

source('minimize.linear_v1.R')

#safe_colorblind_palette
safe4 <- c("#D55E00","#0072B2", "#F0E442", "#CC79A7")

data.subset <- read.table('data.alpha.subset.600_v1.csv',header=T,stringsAsFactors = F,sep=',')

plot(data.subset[,'kurtosis']-3,data.subset[,'skewness']^2,xlim=c(0,6),ylim=c(0,6.5),pch=16)
title('subset of 600 random distributions')

data.alpha <- read.table('data.alpha_v2.csv',header=T,stringsAsFactors = F,sep=',')


# plot for the Box.k.sigma model
x <- c(100,10000,1000000)
y <- c(3.8,6.7,9.4)

q0.25.3sigma <- quantile(data.alpha[,'alpha.plus.3sigma'],0.25)
q0.25.4sigma <- quantile(data.alpha[,'alpha.plus.4sigma'],0.25)
q0.25.5sigma <- quantile(data.alpha[,'alpha.plus.5sigma'],0.25)

q0.50.3sigma <- quantile(data.alpha[,'alpha.plus.3sigma'],0.5)
q0.50.4sigma <- quantile(data.alpha[,'alpha.plus.4sigma'],0.5)
q0.50.5sigma <- quantile(data.alpha[,'alpha.plus.5sigma'],0.5)

q0.75.3sigma <- quantile(data.alpha[,'alpha.plus.3sigma'],0.75)
q0.75.4sigma <- quantile(data.alpha[,'alpha.plus.4sigma'],0.75)
q0.75.5sigma <- quantile(data.alpha[,'alpha.plus.5sigma'],0.75)


x <- c(100,10000,1000000)
y.gauss <- c(1.7239,2.4652,3.207)
y.pearson <- c(3.7975,6.7499,9.37055)
y.exp <- c(4.7497,8.1626,12.4415)

aR.gauss <- minimize.linear(x,y.gauss)
aR.pearson <- minimize.linear(x,y.pearson)
aR.exp <- minimize.linear(x,y.exp)

plot(log(x),c(q0.50.3sigma,q0.50.4sigma,q0.50.5sigma),ylim=c(0,16),xlim=c(4,15),pch=16)
points(log(x),c(q0.25.3sigma,q0.25.4sigma,q0.25.5sigma))
points(log(x),c(q0.75.3sigma,q0.75.4sigma,q0.75.5sigma))
lines(log(x),0.61*log(x)+1)
points(log(x),y.gauss,pch=16)
lines(log(x),0.16*log(x)+1)
points(log(x),y.pearson,pch=16)
lines(log(x),0.81*log(x)+1)
points(log(x),y.exp,pch=16)



coldist <- c('grey','grey','grey','grey','grey','grey')

data.gamma <- data.alpha['gamma'==data.alpha[,1],]
data.inversegamma <- data.alpha['inverse gamma'==data.alpha[,1],]
data.beta <- data.alpha['beta'==data.alpha[,1],]
data.betaprime <- data.alpha['betaprime'==data.alpha[,1],]
data.student <- data.alpha['student'==data.alpha[,1],]
data.pearson <- data.alpha['pearson IV'==data.alpha[,1],]

plot(data.beta[,'kurtosis']-3,data.beta[,'skewness']^2,pch='.',xlim=c(0,6),ylim=c(0,6.5),col=coldist[1])
points(data.betaprime[,'kurtosis']-3,data.betaprime[,'skewness']^2,pch='.',col=coldist[2])
points(data.pearson[,'kurtosis']-3,data.pearson[,'skewness']^2,pch='.',col=coldist[3])
points(data.student[,'kurtosis']-3,data.student[,'skewness']^2,pch='.',col=coldist[4])
points(data.gamma[,'kurtosis']-3,data.gamma[,'skewness']^2,pch='.',col=coldist[5])
points(data.inversegamma[,'kurtosis']-3,data.inversegamma[,'skewness']^2,pch='.',col=coldist[6])




q3s <- quantile(data.alpha[,'alpha.plus.3sigma'])
q3s <- as.numeric(q3s[2:4])

read1 <- data.alpha[,'alpha.plus.3sigma'] < q3s[1]
read2 <- q3s[1] <= data.alpha[,'alpha.plus.3sigma'] & data.alpha[,'alpha.plus.3sigma'] < q3s[2]
read3 <- q3s[2] <= data.alpha[,'alpha.plus.3sigma'] & data.alpha[,'alpha.plus.3sigma'] < q3s[3]
read4 <- q3s[3] <= data.alpha[,'alpha.plus.3sigma']


plot(data.alpha[read1,'kurtosis']-3,data.alpha[read1,'skewness']^2,pch='.',xlim=c(0,6),ylim=c(0,6.5),col=safe4[1])
points(data.alpha[read2,'kurtosis']-3,data.alpha[read2,'skewness']^2,pch='.',col=safe4[2])
points(data.alpha[read3,'kurtosis']-3,data.alpha[read3,'skewness']^2,pch='.',col=safe4[3])
points(data.alpha[read4,'kurtosis']-3,data.alpha[read4,'skewness']^2,pch='.',col=safe4[4])


q4s <- quantile(data.alpha[,'alpha.plus.4sigma'])
q4s <- as.numeric(q4s[2:4])

read1 <- data.alpha[,'alpha.plus.4sigma'] < q4s[1]
read2 <- q4s[1] <= data.alpha[,'alpha.plus.4sigma'] & data.alpha[,'alpha.plus.4sigma'] < q4s[2]
read3 <- q4s[2] <= data.alpha[,'alpha.plus.4sigma'] & data.alpha[,'alpha.plus.4sigma'] < q4s[3]
read4 <- q4s[3] <= data.alpha[,'alpha.plus.4sigma']

plot(data.alpha[read1,'kurtosis']-3,data.alpha[read1,'skewness']^2,pch='.',xlim=c(0,6),ylim=c(0,6.5),col=safe4[1])
points(data.alpha[read2,'kurtosis']-3,data.alpha[read2,'skewness']^2,pch='.',col=safe4[2])
points(data.alpha[read3,'kurtosis']-3,data.alpha[read3,'skewness']^2,pch='.',col=safe4[3])
points(data.alpha[read4,'kurtosis']-3,data.alpha[read4,'skewness']^2,pch='.',col=safe4[4])



q5s <- quantile(data.alpha[,'alpha.plus.5sigma'])
q5s <- as.numeric(q5s[2:4])

read1 <- data.alpha[,'alpha.plus.5sigma'] < q5s[1]
read2 <- q5s[1] <= data.alpha[,'alpha.plus.5sigma'] & data.alpha[,'alpha.plus.5sigma'] < q5s[2]
read3 <- q5s[2] <= data.alpha[,'alpha.plus.5sigma'] & data.alpha[,'alpha.plus.5sigma'] < q5s[3]
read4 <- q5s[3] <= data.alpha[,'alpha.plus.5sigma']


plot(data.alpha[read1,'kurtosis']-3,data.alpha[read1,'skewness']^2,pch='.',xlim=c(0,6),ylim=c(0,6.5),col=safe4[1])
points(data.alpha[read2,'kurtosis']-3,data.alpha[read2,'skewness']^2,pch='.',col=safe4[2])
points(data.alpha[read3,'kurtosis']-3,data.alpha[read3,'skewness']^2,pch='.',col=safe4[3])
points(data.alpha[read4,'kurtosis']-3,data.alpha[read4,'skewness']^2,pch='.',col=safe4[4])



data.0.025 <- read.table('data.0.025_v4.csv',header=T,stringsAsFactors = F,sep=',')
data.0.25 <- read.table('data.0.25_v4.csv',header=T,stringsAsFactors = F,sep=',')
data.0.50 <- read.table('data.0.50_v4.csv',header=T,stringsAsFactors = F,sep=',')
data.0.75 <- read.table('data.0.75_v4.csv',header=T,stringsAsFactors = F,sep=',')
data.0.975 <- read.table('data.0.975_v4.csv',header=T,stringsAsFactors = F,sep=',')

###
###

y.0.025 <- c(data.0.025[1,'Nsize.16384'],data.0.025[2,'Nsize.16384'],data.0.025[3,'Nsize.16384'],data.0.025[4,'Nsize.16384'])
y.25 <- c(data.0.25[1,'Nsize.16384'],data.0.25[2,'Nsize.16384'],data.0.25[3,'Nsize.16384'],data.0.25[4,'Nsize.16384'])
y.50 <- c(data.0.50[1,'Nsize.16384'],data.0.50[2,'Nsize.16384'],data.0.50[3,'Nsize.16384'],data.0.50[4,'Nsize.16384'])
y.75 <- c(data.0.75[1,'Nsize.16384'],data.0.75[2,'Nsize.16384'],data.0.75[3,'Nsize.16384'],data.0.75[4,'Nsize.16384'])
y.0.975 <- c(data.0.975[1,'Nsize.16384'],data.0.975[2,'Nsize.16384'],data.0.975[3,'Nsize.16384'],data.0.975[4,'Nsize.16384'])

y.0.025 <- y.0.025*100
y.25 <- y.25*100
y.50 <- y.50*100
y.75 <- y.75*100
y.0.975 <- y.0.975*100

x <- 1:4
plot(x,y.25,ylim=c(94,100),xlab='Kimber,Hubert,Schwertman,Leys',ylab='percentage captured')
points(x,y.50)
points(x,y.75)
points(x,y.0.025,pch=16)
points(x,y.0.975,pch=16)
title('boxplot (draft for inkscape)')

x <- 1
plot(x,data.0.25[5,'Nsize.16384']*100,ylim=c(99.97,100),xlab='Ritter',ylab='percentage captured')
points(x,data.0.50[5,'Nsize.16384']*100)
points(x,data.0.75[5,'Nsize.16384']*100)
points(x,data.0.025[5,'Nsize.16384']*100,pch=16)
points(x,data.0.975[5,'Nsize.16384']*100,pch=16)
title('boxplot (draft for inkscape)')

Kimber.inf <- data.0.50[1,'Nsize.16384']
Hubert.inf <- data.0.50[2,'Nsize.16384']
Schwertman.inf <- data.0.50[3,'Nsize.16384']
Leys.inf <- data.0.50[4,'Nsize.16384']
Ritter.inf <- data.0.50[5,'Nsize.16384']

data.0.50[1,-1] <- as.numeric(data.0.50[1,-1])-Kimber.inf
data.0.50[2,-1] <- as.numeric(data.0.50[2,-1])-Hubert.inf
data.0.50[3,-1] <- as.numeric(data.0.50[3,-1])-Schwertman.inf
data.0.50[4,-1] <- as.numeric(data.0.50[4,-1])-Leys.inf
data.0.50[5,-1] <- as.numeric(data.0.50[5,-1])-Ritter.inf

data.0.50[,-1] <- round(data.0.50[,-1]*100,digits=2)

Kimber.inf <- round(Kimber.inf*100,digits=2)
Hubert.inf <- round(Hubert.inf*100,digits=2)
Schwertman.inf <- round(Schwertman.inf*100,digits=2)
Leys.inf <- round(Leys.inf*100,digits=2)
Ritter.inf <- round(Ritter.inf*100,digits=2)

# log scale for the sample size
x <- log(c(16,32,64,128,256,512,1024,2048,4096,8192,16384))


plot(x,data.0.50[1,-1],type='l',ylim=c(-6,0),col=safe4[4],lwd=2)
points(x,data.0.50[1,-1],col=safe4[4],pch=16)
lines(x,data.0.50[2,-1],col=safe4[1],lwd=2)
points(x,data.0.50[2,-1],col=safe4[1],pch=16)
lines(x,data.0.50[3,-1],col=safe4[2],lwd=2)
points(x,data.0.50[3,-1],col=safe4[2],pch=16)
lines(x,data.0.50[4,-1],col=safe4[3],lwd=2)
points(x,data.0.50[4,-1],col=safe4[3],pch=16)
lines(x,data.0.50[5,-1],col="black",lwd=2)
points(x,data.0.50[5,-1],col="black",pch=16)
