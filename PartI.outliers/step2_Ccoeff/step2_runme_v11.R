# This script calculates optimum C values (for the boxplot) and beta value (for the MAD rule) on the Pearson family for n=9. The expected percentage of erroneously flagged outliers is 0.1%.

# working folder
path.main <- "~/Science/ctbi11/PartI.outliers/step2_Ccoeff/"
setwd(path.main)

library(gsl) # required for pearsonDS
library(PearsonDS) # for the pearson distribution
library(stats) # for the gamma distribution
library(extraDistr) # for the inverse gamma and betaprime distribution
library(VGAM) # for the weibull, frechet and gumbel distributions

# seed
set.seed(1)

# clarify which package to use
rinvgamma <- extraDistr::rinvgamma
qinvgamma <- extraDistr::qinvgamma
rfrechet <- VGAM::rfrechet
qfrechet <- VGAM::qfrechet
rgumbel <- VGAM::rgumbel
qgumbel <- VGAM::qgumbel

if(0) # estimate the impact of C=34,36,38 on 100 generations of 100 random pearson distributions replicated 1000 times for n=9. Total number of points : 100*100*1000*9 = 9x(10^7). C=36 corresponds to 0.1% of type I error.
{
data.pearson <- read.table('data.alpha_AB_v10.csv',sep=',',header=T)

pearson.family <- c('normal','gamma','inverse gamma','betaprime','pearson IV','student')
data.pearson <- data.pearson[data.pearson[,'distribution'] %in% pearson.family,]

all.size <- 9
all.C <- c(34,36,38)
N.loop <- 1000
N.generations <- 100
N.distributions <- 100
for(k.generation in 1:N.generations)
{
  data.sample <- data.pearson[sample(1:length(data.pearson[,1]),N.distributions,replace=FALSE),]
  
  data.size <- data.frame(generation=rep(k.generation,length(all.size)),sample.size=all.size,stringsAsFactors = F)
  for(k.size in 1:length(all.size))
  {
    N.size <- all.size[k.size]
    for(k.C in 1:length(all.C))
    {
      print(c(k.generation,k.size,k.C))
      
      C <- all.C[k.C]
      
      n.points <- 0
      n.flag <- 0
      for(k.dist in 1:length(data.sample[,1]))
      {
        shape1 <- data.sample[k.dist,'shape1']
        shape2 <- data.sample[k.dist,'shape2']
        A <- data.sample[k.dist,'A']
        B <- data.sample[k.dist,'B']
        char.test <- data.sample[k.dist,'distribution']
        for(k.loop in 1:N.loop)
        {
          if(1)
          {
            if(char.test=='normal')
            {
              d.test <- rnorm(N.size)
            }
            if(char.test=='student')
            {
              d.test <- rt(N.size, df=shape1)
            }
            if(char.test=='gamma')
            {
              d.test <- rgamma(N.size,shape=shape1)
            }
            if(char.test=='inverse gamma')
            {
              d.test <- rinvgamma(N.size,alpha=shape1)
            }
            if(char.test=='pearson IV')
            {
              d.test <- rpearsonIV(N.size,m=shape1,nu=shape2,location=0,scale=1)
            }
            if(char.test=='betaprime')
            {
              d.test <- rbetapr(N.size,shape1=shape1,shape2=shape2)
            }
          }
          
          q0.75 <- quantile(d.test,0.75)
          q0.25 <- quantile(d.test,0.25)
          
          alpha <- B+A*log(N.size)+C/N.size
          
          lower <- q0.25-alpha*(q0.75-q0.25)
          upper <- q0.75+alpha*(q0.75-q0.25)
          
          read_ <- d.test < lower | d.test > upper
          
          n.points <- n.points+N.size
          n.flag <- n.flag+sum(read_)
        }
      }
      
      data.size[k.size,paste0('C',C)] <- n.flag*100/n.points
    }
  }
  
  if(k.generation == 1)
  {
    data.out <- data.size
  }else
  {
    data.out <- rbind(data.out,data.size)
  }
}
write.table(data.out,'Pearson_smallsamples_C34_36_38.csv',row.names = F,sep=',')
}

if(0) # estimate the impact of beta=11,12.5,14 on 100 generations of 100 random pearson distributions replicated 1000 times for n=9. Total number of points : 100*100*1000*9 = 9x(10^7). beta=12.5 corresponds to 0.1% of type I error.
{
  data.pearson <- read.table('data.alpha_AB_v10.csv',sep=',',header=T)
  
  pearson.family <- c('normal','gamma','inverse gamma','betaprime','pearson IV','student')
  data.pearson <- data.pearson[data.pearson[,'distribution'] %in% pearson.family,]
  
  all.size <- 9
  all.beta <- c(11,12.5,14)
  N.loop <- 1000
  N.generations <- 100
  N.distributions <- 100
  for(k.generation in 1:N.generations)
  {
    data.sample <- data.pearson[sample(1:length(data.pearson[,1]),N.distributions,replace=FALSE),]
    
    data.size <- data.frame(generation=rep(k.generation,length(all.size)),sample.size=all.size,stringsAsFactors = F)
    for(k.size in 1:length(all.size))
    {
      N.size <- all.size[k.size]
      for(k.C in 1:length(all.beta))
      {
        print(c(k.generation,k.size,k.C))
        
        beta <- all.beta[k.C]
        
        n.points <- 0
        n.flag <- 0
        for(k.dist in 1:length(data.sample[,1]))
        {
          shape1 <- data.sample[k.dist,'shape1']
          shape2 <- data.sample[k.dist,'shape2']

          char.test <- data.sample[k.dist,'distribution']
          for(k.loop in 1:N.loop)
          {
            if(1)
            {
              if(char.test=='normal')
              {
                d.test <- rnorm(N.size)
              }
              if(char.test=='student')
              {
                d.test <- rt(N.size, df=shape1)
              }
              if(char.test=='gamma')
              {
                d.test <- rgamma(N.size,shape=shape1)
              }
              if(char.test=='inverse gamma')
              {
                d.test <- rinvgamma(N.size,alpha=shape1)
              }
              if(char.test=='pearson IV')
              {
                d.test <- rpearsonIV(N.size,m=shape1,nu=shape2,location=0,scale=1)
              }
              if(char.test=='betaprime')
              {
                d.test <- rbetapr(N.size,shape1=shape1,shape2=shape2)
              }
            }
            
            q0.50 <- quantile(d.test,0.5)

            mad0 <- mad(d.test)
            
            if(mad0 != 0)
            {
            lower <- q0.50-beta*mad0
            upper <- q0.50+beta*mad0
            
            read_ <- d.test < lower | d.test > upper
            
            n.points <- n.points+N.size
            n.flag <- n.flag+sum(read_)
            }
          }
        }
        
        data.size[k.size,paste0('beta',beta)] <- n.flag*100/n.points
      }
    }
    
    if(k.generation == 1)
    {
      data.out <- data.size
    }else
    {
      data.out <- rbind(data.out,data.size)
    }
  }
  write.table(data.out,'Pearson_smallsamples_Beta11_125_14.csv',row.names = F,sep=',')
}

