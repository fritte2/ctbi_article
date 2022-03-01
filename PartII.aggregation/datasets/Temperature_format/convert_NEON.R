setwd("~/Science/PAST/PartII.aggregation_v1/datasets/Temperature_format")

data1 <- read.table('NEON.D17.SJER.DP1.00002.001.000.050.001.SAAT_1min.2020-07.basic.20200803T181022Z.csv',header=T,stringsAsFactors = F,sep=',')
data2 <- read.table('NEON.D17.SJER.DP1.00002.001.000.050.001.SAAT_1min.2020-08.basic.20200902T161325Z.csv',header=T,stringsAsFactors = F,sep=',')

t1 <- as.POSIXct(paste0(substr(data1[,1],1,10),' ',substr(data1[,1],12,19)),tz='UTC')+30
t2 <- as.POSIXct(paste0(substr(data2[,1],1,10),' ',substr(data2[,1],12,19)),tz='UTC')+30


data <- data.frame(time=c(t1,t2),temperature=c(data1[,3],data2[,3]),stringsAsFactors = F)

write.table(data,'temperature_SJER_1min.csv',row.names = F,sep=',')