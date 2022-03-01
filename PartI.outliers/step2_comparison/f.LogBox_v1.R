f.LogBox <- function(d.dist,q0.25,q0.75,k.coeff.outliers)
{
  n <- length(d.dist)
  
  alpha <- k.coeff.outliers*log(n)+1
  
  lower.boundary <- q0.25-alpha*(q0.75-q0.25)
  upper.boundary <- q0.75+alpha*(q0.75-q0.25)
  
  read_ <- lower.boundary <= d.dist & d.dist <= upper.boundary
  
  return(sum(read_))
}