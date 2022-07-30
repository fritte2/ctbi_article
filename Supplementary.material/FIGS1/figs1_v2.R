setwd("~/Science/ctbi12/Supplementary.material/FIGS1/")

data.C <- read.table('Pearson_smallsamples_C34_36_38.csv',header=T,stringsAsFactors = F,sep=',')

data.beta <- read.table('Pearson_smallsamples_Beta11_125_14.csv',header=T,stringsAsFactors = F,sep=',')

k.loop <- 0
for(i in 1:2)
{
  if(i==1)
  {
    data <- data.C
  }
  if(i==2)
  {
    data <- data.beta
  }
  
  for(j in 3:5)
  {
    k.loop <- k.loop+1
    q0.05 <- quantile(data[,j],0.05)
    q0.25 <- quantile(data[,j],0.25)
    q0.50 <- quantile(data[,j],0.5)
    q0.75 <- quantile(data[,j],0.75)
    q0.975 <- quantile(data[,j],0.975)
    
    all.q <- c(q0.05,q0.25,q0.50,q0.75,q0.975)
    
    if(k.loop == 1)
    {
      plot(rep(k.loop,5),all.q,ylim=c(0.02,0.2),xlim=c(0,7),pch='+')
    }else
    {
      points(rep(k.loop,5),all.q,pch='+')
    }
  }
}
abline(h=0.1,col='red')