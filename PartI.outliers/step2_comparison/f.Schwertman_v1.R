f.Schwertman <- function(d.dist,q0.25,q0.50,q0.75,Z,k.n)
{
  lower.boundary <- q0.50-(Z/k.n)*2*(q0.50-q0.25)
  upper.boundary <- q0.50+(Z/k.n)*2*(q0.75-q0.50)
  
  read_ <- lower.boundary <= d.dist & d.dist <= upper.boundary
  
  return(sum(read_))
}