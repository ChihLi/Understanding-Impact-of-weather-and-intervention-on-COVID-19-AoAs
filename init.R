# initialization
init <- function(U){
  
  d <- ncol(U)
  n <- nrow(U)
  
  Theta <- matrix(0, ncol = 2, nrow = n)
  Theta[,1] <- logit.inv(rnorm(n, logit(0.5), 0.05))
  Theta[,2] <- logit.inv(rnorm(n, logit(0.5), 0.05))
  
  ratio <- logit.inv(rnorm(n, logit(0.5), 0.05))
  ratio <- rep(0.5, n)
  
  beta1 <- 0
  beta2 <- 0
  beta.r <- 0
  
  phi1 <- phi2 <- rep(0.5^(1/d), d)
  phi.r <- 0.5
  
  tau <- 1
  tau.r <- 1
  crscor <- 0.9
  
  return(list(Theta = Theta,
              beta1 = beta1, beta2 = beta2,
              phi1 = phi1, phi2 = phi2, 
              tau = tau, tau.r = tau.r,
              crscor = crscor, nug = 1e-4,
              ratio = ratio, beta.r = beta.r, 
              phi.r = phi.r))
}