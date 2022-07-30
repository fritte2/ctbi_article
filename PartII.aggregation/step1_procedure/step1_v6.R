# this scripts plots the aggregation procedure illustrated in FIG. 3 with modified versions of the ctbi scripts.

###
path.main <- "~/Science/ctbi12/PartII.aggregation/step1_procedure/"
file.sources = list.files(pattern="*.R",path=paste0(path.main,'ctbi_showmedian/'),full.names = TRUE)
lapply(file.sources,source)
###

setwd(path.main)
library(data.table)
set.seed(1)

x <- seq(from=as.Date('2017-06-16'),to=as.Date('2022-05-16'),by='1 month')

y <- 5*cos(2*pi*(0:(length(x)-1))/12)+0.9*rnorm(length(x))

xx <- (0:(length(x)-1)/(length(x)-1))*5

y <- y+(0.5*xx^2-4*xx+4)

y[1:24] <- y[1:24]+2

y <- y+12

data <- data.frame(xx=x,y=y)
data[data[,2] < 2,2] <- 2
data[c(3,11,12),2] <- c(7,NA,NA)
data[31,2] <- 17

data <- data[-c(40:45),]

bin.period <- '1 year'
bin.FUN <- 'mean'
bin.center <- NULL
bin.side <- as.Date('2017-06-01')
coeff.outlier <- 'auto'
SCI.min <- Inf
ylim <- c(-Inf,Inf)
bin.max.f.NA <- 0.2

list.main <- ctbi(data,bin.side=bin.side,bin.period=bin.period,bin.center=bin.center,bin.FUN = bin.FUN,bin.max.f.NA = bin.max.f.NA,SCI.min = SCI.min,ylim=ylim,coeff.outlier=coeff.outlier)
data0 <- list.main$data0
data1 <- list.main$data1



