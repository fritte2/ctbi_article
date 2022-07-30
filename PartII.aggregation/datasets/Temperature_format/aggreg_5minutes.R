setwd("~/Science/ctbi12/PartII.aggregation/datasets/Temperature_format/")

library(ctbi)

data <- read.table('temperature_SJER_1min.csv',header=T,stringsAsFactors = F,sep=',')
data[,1] <- as.POSIXct(data[,1])

list.main <- ctbi(data,bin.period='5 min',bin.center=as.POSIXct('2000-01-01 00:02:30'),coeff.outlier = NA)

data1 <- list.main$data1


plot(data[,1],data[,2])
lines(data1[,1],data1[,2],col='red')

write.table(data1,'temperature_SJER.csv',row.names = F,sep=',')
