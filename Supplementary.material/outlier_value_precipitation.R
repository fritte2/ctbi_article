# this script justifies why the outlier level for the precipitation is 1.6 * y.max

path.precipitation <- "~/Science/ctbi_submission/residuals/IN_precipitation/"
setwd(path.precipitation)

library(data.table)
stations <- list.files(pattern='*.csv$')

all.lambda <- c()
for(i in 1:length(stations))
{
  print(i)
  station.i <- fread(stations[i])
  
  y.max <- max(unlist(station.i[,prcp]),na.rm=T)*1.2 # events 20% above the maximum precipitation in the century are considered as impossible
  
  if(nrow(station.i) >= 365*100)
  {
    # randomly select 30 years
  date.1 <- sample(unlist(station.i[,date])[1:(nrow(station.i)-365*30)],1)
  date.2 <- date.1+30*365
  
  
  y.max.30 <- max(unlist(station.i[date.1 <= date & date <= date.2,prcp]),na.rm=T) # calculate the maximum
  
  if(y.max != 0 & y.max.30 != 0)
  {
    all.lambda <- c(all.lambda,y.max/y.max.30)
  }
  }
}

print(paste0('acceptable outlier value for the precipitation : y.outlier = y.max.30years * ',round(mean(all.lambda),digits=1)))


