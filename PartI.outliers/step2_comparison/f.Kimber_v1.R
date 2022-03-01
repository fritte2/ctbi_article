f.Kimber <- function(d.dist,q0.25,q0.50,q0.75)
{
  lower.boundary <- q0.25-3*(q0.50-q0.25)
  upper.boundary <- q0.75+3*(q0.75-q0.50)
  
  read_ <- lower.boundary <= d.dist & d.dist <= upper.boundary
  
  return(sum(read_))
}