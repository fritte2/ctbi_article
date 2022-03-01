setwd("~/Science/hess-2021-609/Supplementary.material/PartII.aggregation")
set.seed(1)

library(ctbi)

safe.colors <- c("#D55E00","#0072B2", "#F0E442", "#CC79A7",'grey','black')

# INPUTS
bin.max.f.NA <- 0.2
bin.period.all <- c('5 months','10 months','15 months','30 months','50 months','100 months')
bin.period.numeric <- c(5,10,15,30,50,100)
bin.center <- NULL
bin.side <- as.Date('2000-01-01')
bin.FUN <- 'mean'
SCI.min <- Inf
k.outliers <- Inf

if(0) # if(1), run a ~4 hour long script to calculate the bias on (1-SS.res/SS.tot). if(0), import the .csv file produced by this script. 
{
number.of.bins <- c(3:10,seq(from=12,to=20,by=2),seq(from=25,to=50,by=5),seq(from=60,to=100,by=10))
number.of.replica <- floor(10^4/number.of.bins)

biased.stacked.cycles.index <- data.frame(number.of.bins=number.of.bins,bin.with.5.points=rep(NA,length(number.of.bins)),bin.with.10.points=rep(NA,length(number.of.bins)),bin.with.15.points=rep(NA,length(number.of.bins)),bin.with.30.points=rep(NA,length(number.of.bins)),bin.with.50.points=rep(NA,length(number.of.bins)),bin.with.100.points=rep(NA,length(number.of.bins)),stringsAsFactors = F)

for(k.period in 1:length(bin.period.all))
{
  bin.period <- bin.period.all[k.period]
  
  for(k.bins in 1:length(number.of.bins))
  {
    all.stacked.cycles.index <- NULL
    for(k.replica in 1:number.of.replica[k.bins])
    {
      print(c(k.bins,length(number.of.bins),k.replica))
      x <- seq(from=as.Date('2000-01-15'),by='1 month',length.out=bin.period.numeric[k.period]*number.of.bins[k.bins])
      data0 <- data.frame(x=x,y=rnorm(bin.period.numeric[k.period]*number.of.bins[k.bins]),stringsAsFactors = F)
      
      list.main <- ctbi(data0,bin.side=bin.side,bin.period=bin.period,bin.center=bin.center,bin.FUN = bin.FUN,bin.max.f.NA = bin.max.f.NA,SCI.min = SCI.min,k.outliers=k.outliers)
      data0 <- list.main$data0
      data1 <- list.main$data1
      stacked.cycles.index <- list.main$SCI+(1/sum(data1[,'index.bin'] > 0))
      
      all.stacked.cycles.index <- c(all.stacked.cycles.index,stacked.cycles.index)
    }
    
    biased.stacked.cycles.index[k.bins,1+k.period] <- mean(all.stacked.cycles.index)
  }
}

write.table(biased.stacked.cycles.index,'biased.stacked.cycles.index.R2_v3.csv',row.names = F,sep=',')
}else
{
  biased.stacked.cycles.index <- read.table('biased.stacked.cycles.index.R2_v2.csv',sep=',',header=T,stringsAsFactors = F)
}

# format of the fit : y = a*(x^b)

# bin with 5 data points
a1 <- lm(log(biased.stacked.cycles.index[,2]) ~ log(biased.stacked.cycles.index[,1]))
b1 <- a1$coefficients[2]
a1 <- a1$coefficients[1]

# bin with 10 data points
a2 <- lm(log(biased.stacked.cycles.index[,3]) ~ log(biased.stacked.cycles.index[,1]))
b2 <- a2$coefficients[2]
a2 <- a2$coefficients[1]

# bin with 15 data points
a3 <- lm(log(biased.stacked.cycles.index[,4]) ~ log(biased.stacked.cycles.index[,1]))
b3 <- a3$coefficients[2]
a3 <- a3$coefficients[1]

# bin with 30 data points
a4 <- lm(log(biased.stacked.cycles.index[,5]) ~ log(biased.stacked.cycles.index[,1]))
b4 <- a4$coefficients[2]
a4 <- a4$coefficients[1]

# bin with 50 data points
a5 <- lm(log(biased.stacked.cycles.index[,6]) ~ log(biased.stacked.cycles.index[,1]))
b5 <- a5$coefficients[2]
a5 <- a5$coefficients[1]

# bin with 100 data points
a6 <- lm(log(biased.stacked.cycles.index[,7]) ~ log(biased.stacked.cycles.index[,1]))
b6 <- a6$coefficients[2]
a6 <- a6$coefficients[1]


a.all <- exp(as.numeric(c(a1,a2,a3,a4,a5,a6)))
b.all <- as.numeric(c(b1,b2,b3,b4,b5,b6))

# a.all is therefore considered as equal to 1 and b.all equal to -1.

plot(biased.stacked.cycles.index[,1],biased.stacked.cycles.index[,2],xlim=c(0,100),pch=16)
points(biased.stacked.cycles.index[,1],biased.stacked.cycles.index[,3],col=safe.colors[1],pch=16)
points(biased.stacked.cycles.index[,1],biased.stacked.cycles.index[,4],col=safe.colors[2],pch=16)
points(biased.stacked.cycles.index[,1],biased.stacked.cycles.index[,5],col=safe.colors[3],pch=16)
points(biased.stacked.cycles.index[,1],biased.stacked.cycles.index[,6],col=safe.colors[4],pch=16)
points(biased.stacked.cycles.index[,1],biased.stacked.cycles.index[,7],col=safe.colors[5],pch=16)
xx <- seq(from=3,to=100,by=1)
lines(xx,1/xx,lwd=1,col='red')




