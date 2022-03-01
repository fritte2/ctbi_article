setwd("~/Science/hess-2021-609/PartII.aggregation/datasets/precipitation_CapeLeuwin")

data <- read.table('PRCP_CapeLeuwin_raw.csv',header=T,stringsAsFactors = F,sep=',')

time <- as.Date(data[,'time'])
PRCP <- as.numeric(data[,'PRCP'])

data <- data.frame(time=time,PRCP=PRCP,stringsAsFactors = F)

# convert occult precipitation to 0
read.low <- data[,2] < 0.01
read.low[is.na(read.low)] <- F
data[as.logical(read.low),2] <- 0

start_ <- as.Date('1920-01-01')
end_ <- as.Date('2019-12-31')

read_ <- start_ <= data[,1] & data[,1] <= end_

data <- data[read_,]

write.table(data,'PRCP_CapeLeuwin.csv',row.names = F,sep=',')