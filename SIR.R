# SIR model
sim.f <- function(x, Theta1, Theta2){
  
  out <- matrix(0,ncol=2,nrow=length(x))

  init <- c(S0, I0, R0)
  N0 <- sum(init)
  names(init) <- c("S", "I", "R")
  
  n <- length(x)
  Iout <- rep(0, n + 1)
  Rout <- rep(0, n + 1)
  Iout[1] <- init["I"]/N0
  Rout[1] <- init["R"]/N0
  for(ii in 2:(n+1)){
    Iout[ii] <- Iout[ii-1] * (1+Theta1[ii-1]-Theta2[ii-1]) - Theta1[ii-1] * Iout[ii-1] * (Iout[ii-1] + Rout[ii-1])
    Rout[ii] <- Rout[ii-1] + Theta2[ii-1] * Iout[ii-1]
  }
  out <- pmax(0,-diff(1-Iout-Rout)* N0)
  
  return(out)
}