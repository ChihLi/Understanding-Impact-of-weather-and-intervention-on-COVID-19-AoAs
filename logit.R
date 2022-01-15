logit <- function(x){
  log(x/(1-x))
}

logit.inv <- function(x){
  exp(x)/(1+exp(x))
}

sqrt.T <- function(crscor){
  A <- matrix(c(1 + sqrt(1-crscor^2), 
                crscor, crscor,
                1 + sqrt(1-crscor^2)), ncol=2, nrow=2)
  A/sqrt(2+2*sqrt(1-crscor^2))
}

compute.jointK <- function(crscor, K1, K2){
  
  A <- sqrt.T(crscor)
  col.n <- ncol(K1)
  row.n <- nrow(K1)
  
  K.big <- matrix(0, ncol=2*col.n, nrow=2*row.n)
  K.big[1:row.n,1:col.n] <- A[1,1]^2*K1+A[1,2]^2*K2
  K.big[(row.n+1):(2*row.n),1:col.n] <- A[1,1]*A[2,1]*K1+A[1,2]*A[2,2]*K2
  K.big[1:row.n, (col.n+1):(2*col.n)] <- A[1,1]*A[1,2]*K1+A[2,1]*A[2,2]*K2
  K.big[(row.n+1):(2*row.n),(col.n+1):(2*col.n)] <- A[2,1]^2*K1+A[2,2]^2*K2   
  
  return(K.big)
}
