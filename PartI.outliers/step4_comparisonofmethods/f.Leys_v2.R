f.Leys <- function(d.dist,q0.50,MAD0)
{
  lower.boundary <- q0.50-3*MAD0
  upper.boundary <- q0.50+3*MAD0
  
  read_ <- lower.boundary <= d.dist & d.dist <= upper.boundary
  
  return(c(sum(!read_),lower.boundary,upper.boundary))
}