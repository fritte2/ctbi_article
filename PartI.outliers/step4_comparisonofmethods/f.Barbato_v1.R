f.Barbato <- function(d.dist,q0.25,q0.75)
{
  n <- length(d.dist)

  alpha <- 1.5*(1+0.1*log(n/10))
  
  lower.boundary <- q0.25-alpha*(q0.75-q0.25)
  upper.boundary <- q0.75+alpha*(q0.75-q0.25)
  
  read_ <- lower.boundary <= d.dist & d.dist <= upper.boundary

  return(c(sum(!read_),lower.boundary,upper.boundary))
}