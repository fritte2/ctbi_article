f.Hubert <- function(d.dist,q0.25,q0.75,MC)
{
  MC.plus <- MC
  if(MC.plus < 0)
  {
    eps.plus <- exp(4*MC.plus)
  }
  else
  {
    eps.plus <- exp(3*MC.plus)
  }
  
  MC.minus <- -MC
  if(MC.minus < 0)
  {
    eps.minus <- exp(4*MC.minus)
  }
  else
  {
    eps.minus <- exp(3*MC.minus)
  }
  
  
  lower.boundary <- q0.25-1.5*eps.minus*(q0.75-q0.25)
  upper.boundary <- q0.75+1.5*eps.plus*(q0.75-q0.25)
  
  read_ <- lower.boundary <= d.dist & d.dist <= upper.boundary
  
  return(c(sum(!read_),lower.boundary,upper.boundary))
}