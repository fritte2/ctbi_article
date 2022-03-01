create.3big.blocks <- function(data0,N.NA.blocks)
{
  N.block.mean <- round(N.NA.blocks/3,digits=0)
  
  i.1 <- sample(10:(round(length(data0[,1])/5)),1)
  
  
  i.2 <- sample((i.1+N.block.mean+10):(length(data0[,1])-2*N.block.mean-20),1)
  
  i.3 <- sample((i.2+N.block.mean+10):(length(data0[,1])-N.block.mean-10),1)
  
  
  read_ <- rep(0,length(data0[,1]))
  
  read_[i.1:(i.1+N.block.mean)] <- 1
  read_[i.2:(i.2+N.block.mean)] <- 1
  read_[i.3:(i.3+N.block.mean)] <- 1
  
  read_ <- as.logical(read_)
  
  index.read <- data0[read_,3]
  return(index.read)
}
