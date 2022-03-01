
setwd("~/Science/ctbi2/hess-2021-609-main_v1.0.0/PartII.aggregation/step1_procedure")
library(ctbi)
set.seed(1)

x <- seq(from=as.Date('2017-06-16'),to=as.Date('2022-05-16'),by='1 month')

y <- 5*cos(2*pi*(0:(length(x)-1))/12)+1.2*rnorm(length(x))

xx <- (0:(length(x)-1)/(length(x)-1))*5

y <- y+(0.5*xx^2-4*xx+4)

y[1:24] <- y[1:24]+2

y <- y+12

data <- data.frame(xx=x,y=y)
data[data[,2] < 2,2] <- 2
data[c(3,11,12),2] <- c(9,NA,NA)
data[31,2] <- 16

data <- data[-c(40:45),]

# plot(x,y)

bin.period <- '1 year'
bin.FUN <- 'mean'
bin.center <- NULL
bin.side <- as.Date('2017-06-01')
k.outliers <- 0.6
SCI.min <- Inf
ylim <- c(-Inf,Inf)
bin.max.f.NA <- 0.2

list.main <- ctbi(data,bin.side=bin.side,bin.period=bin.period,bin.center=bin.center,bin.FUN = bin.FUN,bin.max.f.NA = bin.max.f.NA,SCI.min = SCI.min,ylim=ylim,k.outliers=k.outliers)
data0 <- list.main$data0
data1 <- list.main$data1

plot(data0[,1],data0[,2],pch=16,ylim=c(0,25))
lines(data0[,1],data0[,'long.term'],col='red')
lines(data0[,1],data0[,'long.term']+data0[,'cycle'],col='blue')
abline(v=c(data1[1,'bin.start'],data1[,'bin.end']))

plot(data0[,1],data0[,2]-data0[,'long.term'],pch=16,ylim=c(-6,6))
lines(data0[,1],data0[,'cycle'],col='blue',lty=2)
abline(v=c(data1[1,'bin.start'],data1[,'bin.end']))

plot(data0[,1],data0[,2]-data0[,'long.term']-data0[,'cycle'],pch=16,ylim=c(-6,6))
abline(v=c(data1[1,'bin.start'],data1[,'bin.end']))



