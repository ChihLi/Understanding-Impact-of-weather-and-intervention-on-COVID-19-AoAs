# prior setting
prior <- function(U){
  
  d <- ncol(U)
  
  beta1.mean <- 0
  beta2.mean <- 0
  beta.r.mean <- 0
  
  beta1.var <- 1
  beta2.var <- 1
  beta.r.var <- 1
  
  tau.shape <- 0.01
  tau.rate <- 0.01
  
  tau.r.shape <- 0.01
  tau.r.rate <- 0.01
  
  phi1.b <- phi2.b <- phi.r.b <- 0.1
  cor.b <- 0.1
  
  return(list(beta1.mean = beta1.mean,
              beta2.mean = beta2.mean,
              beta1.var = beta1.var, 
              beta2.var = beta2.var,
              tau.shape = tau.shape,
              tau.rate = tau.rate,
              phi1.b = phi1.b, 
              phi2.b = phi2.b, 
              cor.b = cor.b,
              beta.r.mean = beta.r.mean,
              beta.r.var = beta.r.var, 
              tau.r.shape = tau.r.shape,
              tau.r.rate = tau.r.rate,
              phi.r.b = phi.r.b))
}