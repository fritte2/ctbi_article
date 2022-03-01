create.poisson.blocks <- function(data0,N,N.blocks)
{
  # block are considered here as NA but they can be outliers as well

  epsilon.blocks <- -1+length(data0[,1])/N.blocks

  intervals <- rpois(N.blocks,epsilon.blocks)
  
  # the first data points contains data.
  read_ <- rep(0,N)

  k <- 1
  for(i in 1:length(intervals))
  {
      k <- k+1 # jump to the next point
      k <- k+intervals[i]
      if(k < length(data0[,1]))
      {
      read_[k] <- 1
      }
  }
  
  read_ <- as.logical(read_)
  
  index.read <- data0[read_,3]
  return(index.read)
}
