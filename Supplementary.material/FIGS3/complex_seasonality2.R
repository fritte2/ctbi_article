rm(list = ls())

library(devtools)
devtools::install_github("bpbond/cosore")
library(cosore)
library(data.table)
library(ctbi)
data.all <- setDT(csr_dataset('d20190919_MIGLIAVACCA')$data)

# calculate the mean time of CSR_TIMESTAMP_BEGIN and CSR_TIMESTAMP_END, and then remove them.
CSR.mean <- do.call(c, Map(function(x, y) mean(c(x, y)),data.all[,CSR_TIMESTAMP_BEGIN],data.all[,CSR_TIMESTAMP_END]))

# store time zone
tzone.CSR <- copy(attr(CSR.mean,'tzone'))

# convert to UTC
attr(CSR.mean,'tzone') <- 'UTC'
data.all[,':='(CSR_TIME = CSR.mean,CSR_TIMESTAMP_BEGIN = NULL,CSR_TIMESTAMP_END = NULL )]

# There are 8 different ports (CSR). Each port performs a measurement every 32 minutes, and there is a 4 min gap between two sensors.
# For one port, a day is therefore 45 measurements of 32 minutes.
# Protocol: aggregate every day (45 points of measurements), then every month (28 to 31 points).

bin.side <- as.POSIXct('2000-01-01 00:00:00',tz='UTC')
bin.period.successive <- c('1 day','1 month')
coeff.outlier.successive <- c(NA,NA) # change to c('auto','auto') to see the impact on the aggregation
bin.max.f.NA.successive <- c(0.2,0.3) # days with less than 20% of data are rejected, and months with less than 70% of data are rejected.
SCI.min <- NA # do not impute data
bin.FUN <- 'sum'

data.agg.temp <- copy(data.all)
for(k.agg in 1:2) # perform 2 successive aggregations
{
  bin.period <- bin.period.successive[k.agg]
  coeff.outlier <- coeff.outlier.successive[k.agg]
  bin.max.f.NA <- bin.max.f.NA.successive[k.agg]

  data.agg <- copy(data.agg.temp)
  data.agg <- data.agg[order(CSR_TIME),]
  rm(data.agg.temp)

for(k.port in 1:8) # read the 8 ports
{
# STEP 1, ignore outliers and separate the data into two parts : low variability around the median, and high variability around the median.
list.ctbi <- ctbi(data.agg[CSR_PORT == k.port,c('CSR_TIME','CSR_FLUX_CO2'),with=F],bin.FUN='median',SCI.min=SCI.min,bin.side=bin.side,bin.period=bin.period,coeff.outlier = NA,bin.max.f.NA=bin.max.f.NA)
data0 <- list.ctbi$data0
data1 <- list.ctbi$data1

index.bins.low <- unlist(data1[mad.CSR_FLUX_CO2 < median(mad.CSR_FLUX_CO2,na.rm=T),index.bin],use.names = F)
index.bins.high <- unlist(data1[mad.CSR_FLUX_CO2 >= median(mad.CSR_FLUX_CO2,na.rm=T),index.bin],use.names = F)
rm(data1)

data.low <- data0[index.bin %in% index.bins.low,c('CSR_TIME','CSR_FLUX_CO2'),with=F]
data.high <- data0[index.bin %in% index.bins.high,c('CSR_TIME','CSR_FLUX_CO2'),with=F]
rm(data0)

# STEP 2, aggregate separately the two datasets
list.ctbi <- ctbi(data.low[,c('CSR_TIME','CSR_FLUX_CO2'),with=F],bin.FUN=bin.FUN,SCI.min=SCI.min,bin.side=bin.side,bin.period=bin.period,coeff.outlier = coeff.outlier,bin.max.f.NA=bin.max.f.NA)
data0.low <- list.ctbi$data0
data1.low <- list.ctbi$data1

# uncomment to visualize
# ctbi.plot(list.ctbi)

list.ctbi <- ctbi(data.high[,c('CSR_TIME','CSR_FLUX_CO2'),with=F],bin.FUN=bin.FUN,SCI.min=SCI.min,bin.side=bin.side,bin.period=bin.period,coeff.outlier = coeff.outlier,bin.max.f.NA=bin.max.f.NA)
data0.high <- list.ctbi$data0
data1.high <- list.ctbi$data1

# uncomment to visualize
#ctbi.plot(list.ctbi)

# step 3, merge datasets
data0.temp <- rbind(data0.low[,c('CSR_TIME','CSR_FLUX_CO2'),with=F],data0.high[,c('CSR_TIME','CSR_FLUX_CO2'),with=F])
data0.temp[,CSR_PORT := k.port]

data1.temp <- rbind(data1.low[index.bin > 0,c('CSR_TIME','CSR_FLUX_CO2'),with=F],data1.high[index.bin > 0,c('CSR_TIME','CSR_FLUX_CO2'),with=F])
data1.temp[,CSR_PORT := k.port]

if(k.port == 1)
{
  data.agg.temp <- copy(data1.temp)
}else
{
  data.agg.temp <- rbind(data.agg.temp,data1.temp)
}

if(k.agg == 1)
{
  if(k.port == 1)
  {
    data.clean.32.min <- copy(data0.temp)
  }else
  {
    data.clean.32.min <- rbind(data.clean.32.min,data0.temp)
  }
}

if(k.agg == 2)
{
  if(k.port == 1)
  {
    data.clean.1day <- copy(data0.temp)
    data.clean.1month <- copy(data1.temp)
  }else
  {
    data.clean.1day <- rbind(data.clean.1day,data0.temp)
    data.clean.1month <- rbind(data.clean.1month,data1.temp)
  }
}

per.low <- round(sum(data1.low[n.points != 0,'n.outliers'])*100/sum(data1.low[n.points != 0,'n.points']),digits=1)
per.high <- round(sum(data1.high[n.points != 0,'n.outliers'])*100/sum(data1.high[n.points != 0,'n.points']),digits=1)

if(k.port == 1)
{
  print('###################################################')
}
print(paste0('######### Aggregation ',bin.period,' port ',k.port,' #########'))
print('Percentage of outliers removed:')
print(paste0('Low variability : ',per.low,' %'))
print(paste0('High variability : ',per.high,' %'))
}
}

col_ <- c('black','red','blue','green','cyan','orange','purple','yellow')
for(i in 1:8)
{
  data.temp <- data.clean.1day[CSR_PORT == i,]
  data.temp <- data.temp[order(CSR_TIME),]
  if(i == 1)
  {
    plot(data.temp[,CSR_TIME],data.temp[,CSR_FLUX_CO2],ylim=c(-400,400),type='l',col=col_[i],xlab='time (daily)',ylab='CO2 Flux (sum)',main='CO2 flux (daily sum) of the 8 ports')
  }else
  {
    lines(data.temp[,CSR_TIME],data.temp[,CSR_FLUX_CO2],col=col_[i])
  }
}


