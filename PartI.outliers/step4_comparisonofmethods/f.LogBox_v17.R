f.LogBox <- function(d.dist,coeff.outlier,q0.125,q0.25,q0.375,q0.625,q0.75,q0.875)
{
  n <- length(d.dist)
  
  if(sum(coeff.outlier == 'auto') == 1) # A,B,C are not provided, we calculate them
  {
    m.plus <- (q0.875-q0.625)/(q0.75-q0.25)
    m.minus <- (q0.375-q0.125)/(q0.75-q0.25)
    
    m.star <- round(max(c(m.minus,m.plus))-0.6165,digits=4)
    
    if(m.star < 0){m.star <- 0}
    if(m.star > 2){m.star <- 2}

      # version 13
      a1 <- 0.2294
      a2 <- 2.9416
      a3 <- -0.0512
      a4 <- -0.0684
      A <- round(a1*exp(a2*m.star+a3*m.star^2+a4*m.star^3),digits=2)
      
      # version 12
      b1 <- 1.0585
      b2 <- 15.6960
      b3 <- -17.3618
      b4 <- 28.3511
      b5 <- -11.4726
      B <- round(b1+b2*m.star+b3*(m.star^2)+b4*(m.star^3)+b5*(m.star^4),digits=2)
    
  }else
  {
    A <- coeff.outlier[1]
    B <- coeff.outlier[2]
  }
  
  C <- 36
  alpha <- A*log(n)+B+C/n
  
  lower.boundary <- q0.25-alpha*(q0.75-q0.25)
  upper.boundary <- q0.75+alpha*(q0.75-q0.25)
  
  read_ <- lower.boundary <= d.dist & d.dist <= upper.boundary
  
  coeff.outlier <- c(A,B)
  
  return(c(sum(!read_),lower.boundary,upper.boundary,coeff.outlier))
}