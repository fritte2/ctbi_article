minimize.linear <- function(x,y)
{

  N <- 10^4
  min_ <- 0.1
  max_ <- 1
  a <- runif(N,min=min_,max=max_)
  
  SS.tot <- sum((y-mean(y))^2)
  
  eps <- Inf
  
  for(i in 1:N)
  {
    SS.res <- sum((y-(a[i]*log(x)+1))^2)
    
    ratio <- SS.res/SS.tot
    
    if(ratio < eps)
    {
      a.final <- a[i]
      eps <- ratio
      R.squared <- round(1-ratio,digits=3)
    }
  }
  
  a.final <- round(a.final,digits=2)
  
  return(c(a.final,R.squared))
}