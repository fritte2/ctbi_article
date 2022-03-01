setwd("~/Science/hess-2021-609/Supplementary.material/PartI.outliers")

library(extraDistr)
library(gld)

safe6 <- c("#D55E00","#0072B2", "#F0E442", "#CC79A7",'grey','black')

data.GLD <- read.table('data.gld_v3.csv',sep=',',header=T,stringsAsFactors = F)
data.600 <- read.table('data.alpha.subset.600_v1.csv',sep=',',header=T,stringsAsFactors = F)

x <- 0:6
plot(data.GLD[,'kurtosis']-3,data.GLD[,'skewness']^2,pch='.',col='black',xlim=c(0,6),ylim=c(0,6.5))
lines(x,(4/5)*(x+2))
title('distributions from GLD')


read.beta <- 'beta' == data.600[,1]
read.betaprime <- 'betaprime' == data.600[,1]
read.gamma <- 'gamma' == data.600[,1]
read.inversegamma <- 'inverse gamma' == data.600[,1]
read.pearson <- 'pearson IV' == data.600[,1]
read.student <- 'student' == data.600[,1]

x <- 0:6
plot(data.600[read.beta,'kurtosis']-3,data.600[read.beta,'skewness']^2,pch=16,col=safe6[1],xlim=c(0,6),ylim=c(0,6.5))
points(data.600[read.betaprime,'kurtosis']-3,data.600[read.betaprime,'skewness']^2,pch=16,col=safe6[2])
points(data.600[read.gamma,'kurtosis']-3,data.600[read.gamma,'skewness']^2,pch=16,col=safe6[3])
points(data.600[read.inversegamma,'kurtosis']-3,data.600[read.inversegamma,'skewness']^2,pch=16,col=safe6[4])
points(data.600[read.pearson,'kurtosis']-3,data.600[read.pearson,'skewness']^2,pch=16,col=safe6[5])
points(data.600[read.student,'kurtosis']-3,data.600[read.student,'skewness']^2,pch=16,col=safe6[6])
lines(x,(4/5)*(x+2))
title('distributions from 600')


alpha.GLD.3sigma <- round(median(data.GLD[,'alpha.plus.3sigma']),digits=1)
alpha.GLD.4sigma <- round(median(data.GLD[,'alpha.plus.4sigma']),digits=1)
alpha.GLD.5sigma <- round(median(data.GLD[,'alpha.plus.5sigma']),digits=1)


# betaprime, shape1 and shape2 close to S^2 = 3, k.ex = 5

case.600 <- data.600[427,]
shape1.betaprime <- case.600[1,'shape1']
shape2.betaprime <- case.600[1,'shape2']
a.m.3sigma.betaprime <- case.600[1,'alpha.minus.3sigma']
a.m.4sigma.betaprime <- case.600[1,'alpha.minus.4sigma']
a.m.5sigma.betaprime <- case.600[1,'alpha.minus.5sigma']
a.p.3sigma.betaprime <- case.600[1,'alpha.plus.3sigma']
a.p.4sigma.betaprime <- case.600[1,'alpha.plus.4sigma']
a.p.5sigma.betaprime <- case.600[1,'alpha.plus.5sigma']

# GLD, lambda3 and lambda4 close to S^2 = 3, k.ex = 5
case.GLD <- data.GLD[2446,]
lambda3.GLD <- case.GLD[1,'lambda3']
lambda4.GLD <- case.GLD[1,'lambda4']
a.m.3sigma.GLD <- case.GLD[1,'alpha.minus.3sigma']
a.m.4sigma.GLD <- case.GLD[1,'alpha.minus.4sigma']
a.m.5sigma.GLD <- case.GLD[1,'alpha.minus.5sigma']
a.p.3sigma.GLD <- case.GLD[1,'alpha.plus.3sigma']
a.p.4sigma.GLD <- case.GLD[1,'alpha.plus.4sigma']
a.p.5sigma.GLD <- case.GLD[1,'alpha.plus.5sigma']

n.size <- 10^6

dist.betaprime <- rbetapr(n.size,shape1=shape1.betaprime,shape2=shape2.betaprime)
dist.GLD <- rgl(n.size,lambda2 = 1,lambda3 = lambda3.GLD, lambda4 = lambda4.GLD)

d1 <- (dist.betaprime-mean(dist.betaprime))/sd(dist.betaprime)
d2 <- (dist.GLD-mean(dist.GLD))/sd(dist.GLD)

kurto1 <- mean(d1^4)
kurto2 <- mean(d2^4)
skew1 <- mean(d1^3)
skew2 <- mean(d2^3)




l1.3sigma <- quantile(d1,0.25)-a.m.3sigma.betaprime*(quantile(d1,0.75)-quantile(d1,0.25))
l1.4sigma <- quantile(d1,0.25)-a.m.4sigma.betaprime*(quantile(d1,0.75)-quantile(d1,0.25))
l1.5sigma <- quantile(d1,0.25)-a.m.5sigma.betaprime*(quantile(d1,0.75)-quantile(d1,0.25))

u1.3sigma <- quantile(d1,0.75)+a.p.3sigma.betaprime*(quantile(d1,0.75)-quantile(d1,0.25))
u1.4sigma <- quantile(d1,0.75)+a.p.4sigma.betaprime*(quantile(d1,0.75)-quantile(d1,0.25))
u1.5sigma <- quantile(d1,0.75)+a.p.5sigma.betaprime*(quantile(d1,0.75)-quantile(d1,0.25))

u1.test1 <- quantile(d1,0.75)+6.7*(quantile(d1,0.75)-quantile(d1,0.25))
u1.test2 <- quantile(d1,0.75)+9.4*(quantile(d1,0.75)-quantile(d1,0.25))

l2.3sigma <- quantile(d2,0.25)-a.m.3sigma.GLD*(quantile(d2,0.75)-quantile(d2,0.25))
l2.4sigma <- quantile(d2,0.25)-a.m.4sigma.GLD*(quantile(d2,0.75)-quantile(d2,0.25))
l2.5sigma <- quantile(d2,0.25)-a.m.5sigma.GLD*(quantile(d2,0.75)-quantile(d2,0.25))

u2.3sigma <- quantile(d2,0.75)+a.p.3sigma.GLD*(quantile(d2,0.75)-quantile(d2,0.25))
u2.4sigma <- quantile(d2,0.75)+a.p.4sigma.GLD*(quantile(d2,0.75)-quantile(d2,0.25))
u2.5sigma <- quantile(d2,0.75)+a.p.5sigma.GLD*(quantile(d2,0.75)-quantile(d2,0.25))

hist(d1,xlim=c(-5,16),100,ylim=c(0,1.5*10^5))
abline(v=c(l1.3sigma,u1.3sigma))
abline(v=c(l1.4sigma,u1.4sigma))
abline(v=c(l1.5sigma,u1.5sigma))
# abline(v=c(u1.test1,u1.test2),col='red',lwd=3)

hist(d2,xlim=c(-5,16),30,ylim=c(0,6*10^5))
abline(v=c(l2.3sigma,u2.3sigma))
abline(v=c(l2.4sigma,u2.4sigma))
abline(v=c(l2.5sigma,u2.5sigma))


