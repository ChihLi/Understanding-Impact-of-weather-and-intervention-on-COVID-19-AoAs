library(plgp)
library(MASS)
# initialization
init <- function(U){
  
  d <- ncol(U)
  n <- nrow(U)
  
  Theta <- matrix(0, ncol = 2, nrow = n)
  Theta[,1] <- logit.inv(rnorm(n, logit(0.5), 0.05))
  Theta[,2] <- logit.inv(rnorm(n, logit(0.5), 0.05))
  
  ratio <- logit.inv(rnorm(n, logit(0.5), 0.05))
  
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

# prior
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


set.seed(1)


# setup
nsamples <- 2000
nburnin <- 2000

init.ls <- init(U)
Theta <- init.ls$Theta
psi.ls <- init.ls$psi.ls
g <- init.ls$g
prior.ls <- prior(U)
x.max <- length(x)

acceptCount <- acceptCount.ratio <- 0  
acceptCountPhi1 <- acceptCountPhi2 <- rep(0, ncol(U))
acceptCountPhi.r <- 0
acceptCountRho <- 0

cTheta <- cRatio <- 0.05
cPhi1 <- cPhi2 <- rep(0.5, ncol(U))
cPhi.r <- 0.5
cRho <- 0.5

# priors
beta1.prior.mean <- prior.ls$beta1.mean
beta2.prior.mean <- prior.ls$beta2.mean
beta1.prior.var <- prior.ls$beta1.var
beta2.prior.var <- prior.ls$beta2.var
beta.r.prior.mean <- prior.ls$beta.r.mean
beta.r.prior.var <- prior.ls$beta.r.var

tau.prior.shape <- prior.ls$tau.shape
tau.prior.rate <- prior.ls$tau.rate
tau.r.prior.shape <- prior.ls$tau.r.shape
tau.r.prior.rate <- prior.ls$tau.r.rate

phi1.b <- prior.ls$phi1.b
phi2.b <- prior.ls$phi2.b
phi.r.b <- prior.ls$phi.r.b
cor.b <- prior.ls$cor.b

# dimension
n <- length(y)
d <- ncol(U)

# initial values
Theta <- init.ls$Theta
ratio <- init.ls$ratio
beta1 <- init.ls$beta1
beta2 <- init.ls$beta2
beta.r <- init.ls$beta.r
phi1 <- init.ls$phi1
phi2 <- init.ls$phi2
phi.r <- init.ls$phi.r
tau <- init.ls$tau
tau.r <- init.ls$tau.r
nug <- init.ls$nug
crscor <- init.ls$crscor

# store MC samples
Theta.sample <- array(0, c(n, 2, nsamples/2))
ratio.sample <- matrix(0, ncol=n, nrow=nsamples/2)
tau.sample <- tau.r.sample <- rep(0, nsamples/2)
phi1.sample <- phi2.sample <- matrix(0, ncol = d, nrow = nsamples/2)
phi.r.sample <- rep(0, nrow = nsamples/2)
crscor.sample <- rep(0, nsamples/2)
beta1.sample <- beta2.sample <- beta.r.sample <- rep(0, nsamples/2)
testTheta.sample <- array(0, c(n.test, 2, nsamples/2))
testRatio.sample <- matrix(0, ncol = n.test, nrow = nsamples/2)
testY.sample <- matrix(0, ncol = n.test, nrow = nsamples/2)
X <- matrix(1, nrow = n)
K1 <- covar.sep(U, d=-1/(4*log(phi1)), g=nug)
K2 <- covar.sep(U, d=-1/(4*log(phi2)), g=nug)
mu1 <- X %*% beta1
mu2 <- X %*% beta2
K.r <- covar.sep(matrix(x/x.max,ncol=1), d=-1/(4*log(phi.r)), g=nug)
mu.r <- X %*% beta.r

# build covariance matrix
K.big <- compute.jointK(crscor, K1, K2)

# covariance inverse
Kchol <- chol(K.big)
Kchol.inv <- solve(Kchol)
K.inv <- Kchol.inv %*% t(Kchol.inv)

Kchol.r <- chol(K.r)
Kchol.r.inv <- solve(Kchol.r)
K.r.inv <- Kchol.r.inv %*% t(Kchol.r.inv)

for(i in 1:(nburnin + nsamples)){

  om_temp <- logit(Theta)
  temp <- logit.inv(om_temp)
  
  om.r_temp <- logit(ratio)
  temp.r <- logit.inv(om.r_temp)
  
  for(ii in 1:5){
    
    lambdai <- sim.f(x, temp[,1], temp[,2])
    
    om.r_star <- cRatio * sqrt(tau.r) * t(Kchol.r) %*% rnorm(n) + om.r_temp #* sqrt(1-cTheta^2)
    
    ratio.old <- logit.inv(om.r_temp)
    ratio.star <- logit.inv(om.r_star)
    
    if(any(ratio.star > 0.999| ratio.star < 0.001)) next  # if the thetas are not reasonable, then reject
    
    logRatioStar <- sum(dpois(y, ratio.star*lambdai, log = TRUE)) -
      0.5/tau.r * t(om.r_star - mu.r) %*% K.r.inv %*% (om.r_star - mu.r)
    
    logRatio <- sum(dpois(y, ratio.old*lambdai, log = TRUE)) -
      0.5/tau.r * t(om.r_temp - mu.r) %*% K.r.inv %*% (om.r_temp - mu.r)
    
    logR <- logRatioStar - logRatio
    
    if(exp(logR) >= 1){
      om.r_temp <- om.r_star
      temp.r <- logit.inv(om.r_temp)
      acceptCount.ratio <- acceptCount.ratio + 1
    }else if(runif(1) <= exp(logR)){
      om.r_temp <- om.r_star
      temp.r <- logit.inv(om.r_temp)
      acceptCount.ratio <- acceptCount.ratio + 1
    }
  }
  om.r <- om.r_temp
  ratio <- logit.inv(om.r)
  
  # MH sampler: update Theta
  for(ii in 1:5){
    om_star <- cTheta * sqrt(tau) * t(Kchol) %*% rnorm(2*n) + c(om_temp[,1], om_temp[,2]) 
    om_star <- cbind(om_star[1:n], om_star[(n+1):(2*n)])
    
    Theta.star <- logit.inv(om_star)
    if(any(Theta.star > 0.999| Theta.star < 0.001)) next  # if the thetas are not reasonable, then reject
    
    lambdai.old <- sim.f(x, temp[,1], temp[,2])
    lambdai.star <- sim.f(x, Theta.star[,1], Theta.star[,2])
    
    logThetaStar <- sum(dpois(y, ratio*lambdai.star, log = TRUE)) -
      0.5/tau * t(c(om_star[,1], om_star[,2]) - c(mu1, mu2)) %*% K.inv %*% (c(om_star[,1], om_star[,2]) - c(mu1, mu2)) 
    logTheta <- sum(dpois(y, ratio*lambdai.old, log = TRUE)) -
      0.5/tau * t(c(om_temp[,1], om_temp[,2]) - c(mu1, mu2)) %*% K.inv %*% (c(om_temp[,1], om_temp[,2]) - c(mu1, mu2)) 
    
    logR <- logThetaStar - logTheta
    
    if(exp(logR) >= 1){
      om_temp <- om_star
      temp <- logit.inv(om_temp)
      acceptCount <- acceptCount + 1
    }else if(runif(1) <= exp(logR)){
      om_temp <- om_star
      temp <- logit.inv(om_temp)
      acceptCount <- acceptCount + 1
    }
  }
  om <- om_temp
  Theta <- logit.inv(om)
  
  # Draw beta.r 
  X <- matrix(1, nrow = n)
  beta.post.mean <-
    solve((t(X) %*% K.r.inv %*% X)/tau.r + diag(1/beta.r.prior.var)) %*%
    (beta.r.prior.mean/beta.r.prior.var + (t(X) %*% K.r.inv %*% logit(ratio))/tau.r)
  beta.post.var <- solve((t(X) %*% K.r.inv %*% X)/tau.r + diag(1/beta.r.prior.var))
  beta.r <- t(chol(beta.post.var)) %*% rnorm(1) + beta.post.mean
  beta.r <- drop(beta.r)
  
  mu.r <- X %*% beta.r
  
  # Draw tau.r
  tau.r.post.shape <- tau.r.prior.shape + n/2
  tau.r.post.rate <- tau.r.prior.rate + c(1/2* t(logit(ratio) - mu.r) %*% K.r.inv %*% (logit(ratio) - mu.r))
  
  tau.r <- 1/rgamma(1, tau.r.post.shape, tau.r.post.rate)
  
  # Draw phi.r
  nu.r <- log(-log(phi.r))
  nu.r.star <- nu.r
  nu.r.star <- nu.r.star + cPhi.r * rnorm(1)
  phi.r.star <- exp(-exp(nu.r.star))
  if(phi.r.star > 0.1 & phi.r.star < 0.9){
    K.r.star <- covar.sep(matrix(x/x.max,ncol=1), d=-1/(4*log(phi.r.star)), g=nug)
    KstarKt <- chol(K.r.star)
    KstarKt.inv <- solve(KstarKt)
    K.r.star.inv <- KstarKt.inv %*% t(KstarKt.inv)
    detK.r.starInv_apx <- -2*sum(log(diag(KstarKt)))
    
    detK.r.Inv_apx <- -2*sum(log(diag(Kchol.r)))
    
    log_rhotStar <- 0.5*detK.r.starInv_apx - 
      0.5/tau.r * c( t(om.r - mu.r) %*% K.r.star.inv %*% (om.r - mu.r)) +
      (phi.r.b - 1)* sum(log(1 - phi.r.star)) + sum(nu.r.star - exp(nu.r.star))
    
    log_rhot <- 0.5*detK.r.Inv_apx - 
      0.5/tau.r * c( t(om.r - mu.r) %*% K.r.inv %*% (om.r - mu.r)) +
      (phi.r.b - 1)*sum(log(1 - phi.r)) + sum(nu.r - exp(nu.r))
    
    # Acceptance crieterion
    logR <- log_rhotStar - log_rhot
    xcptProb <- exp(logR)
    
    if(runif(1) <= xcptProb) {
      phi.r <- phi.r.star
      K.r <- K.r.star
      Kchol.r <- KstarKt
      Kchol.r.inv <- KstarKt.inv
      K.r.inv <- K.r.star.inv
      acceptCountPhi.r <- acceptCountPhi.r + 1
    }
  } 
  
  
  
  # Draw beta1 and beta2
  X <- matrix(1, nrow = n)
  Theta.vt <- c(Theta[,1], Theta[,2])
  X.big <- cbind(X,matrix(0,nrow=n,ncol=1))
  X.big <- rbind(X.big, cbind(matrix(0,nrow=n,ncol=1),X))
  beta.prior.var <- c(beta1.prior.var, beta2.prior.var)
  beta.prior.mean <- c(beta1.prior.mean, beta2.prior.mean)
  beta.post.mean <-
    solve((t(X.big) %*% K.inv %*% X.big)/tau + diag(1/beta.prior.var)) %*%
    (beta.prior.mean/beta.prior.var + (t(X.big) %*% K.inv %*% logit(Theta.vt))/tau)
  beta.post.var <- solve((t(X.big) %*% K.inv %*% X.big)/tau + diag(1/beta.prior.var))
  beta <- t(chol(beta.post.var)) %*% rnorm(2) + beta.post.mean
  beta1 <- beta[1]
  beta2 <- beta[2]
  
  mu1 <- X %*% beta1
  mu2 <- X %*% beta2
  
  tau.post.shape <- tau.prior.shape + n/2*2
  tau.post.rate <- tau.prior.rate + c(1/2* t(logit(c(Theta[,1],Theta[,2])) - c(mu1, mu2)) %*% K.inv %*% (logit(c(Theta[,1],Theta[,2])) - c(mu1, mu2)))
  tau <- 1/rgamma(1, tau.post.shape, tau.post.rate)
  
  for(jj in 1:length(phi1)){
    nu <- log(-log(phi1))
    nu.star <- nu
    nu.star[jj] <- nu.star[jj] + cPhi1[jj] * rnorm(1)
    phi1.star <- exp(-exp(nu.star))
    if(phi1.star[jj] < 0.1^(1/d) | phi1.star[jj] > 0.9^(1/d)) next
    
    K1.star <- covar.sep(U, d=-1/(4*log(phi1.star)), g=nug)
    K.big <- compute.jointK(crscor, K1.star, K2)
    KstarKt <- chol(K.big)
    KstarKt.inv <- solve(KstarKt)
    Kstar.inv <- KstarKt.inv %*% t(KstarKt.inv)
    detKstarInv_apx <- -2*sum(log(diag(KstarKt)))
    
    detKInv_apx <- -2*sum(log(diag(Kchol)))
    
    log_rhotStar <- 0.5*detKstarInv_apx -
      0.5/tau * c( t(c(om[,1], om[,2]) - c(mu1, mu2)) %*% Kstar.inv %*% (c(om[,1], om[,2]) - c(mu1, mu2))) +
      (phi1.b - 1)* sum(log(1 - phi1.star)) + sum(nu.star - exp(nu.star))
    
    log_rhot <- 0.5*detKInv_apx -
      0.5/tau * c( t(c(om[,1], om[,2]) - c(mu1, mu2)) %*% K.inv %*% (c(om[,1], om[,2]) - c(mu1, mu2))) +
      (phi1.b - 1)*sum(log(1 - phi1)) + sum(nu - exp(nu))
    
    # Acceptance crieterion
    logR <- log_rhotStar - log_rhot
    xcptProb <- exp(logR)
    
    if(runif(1) <= xcptProb) {
      phi1 <- phi1.star
      K1 <- K1.star
      Kchol <- KstarKt
      Kchol.inv <- KstarKt.inv
      K.inv <- Kstar.inv
      acceptCountPhi1[jj] <- acceptCountPhi1[jj] + 1
    }
  }
  
  for(jj in 1:length(phi2)){
    nu <- log(-log(phi2))
    nu.star <- nu
    nu.star[jj] <- nu.star[jj] + cPhi2[jj] * rnorm(1)
    phi2.star <- exp(-exp(nu.star))
    if(phi2.star[jj] < 0.1^(1/d) | phi2.star[jj] > 0.9^(1/d)) next
    
    K2.star <- covar.sep(U, d=-1/(4*log(phi2.star)), g=nug)
    K.big <- compute.jointK(crscor, K1, K2.star)
    KstarKt <- chol(K.big)
    KstarKt.inv <- solve(KstarKt)
    Kstar.inv <- KstarKt.inv %*% t(KstarKt.inv)
    detKstarInv_apx <- -2*sum(log(diag(KstarKt)))
    
    detKInv_apx <- -2*sum(log(diag(Kchol)))
    
    log_rhotStar <- 0.5*detKstarInv_apx -
      0.5/tau * c( t(c(om[,1], om[,2]) - c(mu1, mu2)) %*% Kstar.inv %*% (c(om[,1], om[,2]) - c(mu1, mu2))) +
      (phi2.b - 1)* sum(log(1 - phi2.star)) + sum(nu.star - exp(nu.star))
    
    log_rhot <- 0.5*detKInv_apx -
      0.5/tau * c( t(c(om[,1], om[,2]) - c(mu1, mu2)) %*% K.inv %*% (c(om[,1], om[,2]) - c(mu1, mu2))) +
      (phi2.b - 1)*sum(log(1 - phi2)) + sum(nu - exp(nu))
    
    # Acceptance crieterion
    logR <- log_rhotStar - log_rhot
    xcptProb <- exp(logR)
    
    if(runif(1) <= xcptProb) {
      phi2 <- phi2.star
      K2 <- K2.star
      Kchol <- KstarKt
      Kchol.inv <- KstarKt.inv
      K.inv <- Kstar.inv
      acceptCountPhi2[jj] <- acceptCountPhi2[jj] + 1
    }
  }
  
  
  # MH sampler: update crscor
  nu <- log(-log(crscor))
  nu.star <- nu + cRho * rnorm(length(nu))
  crscor.star <- exp(-exp(nu.star))
  
  K.big <- compute.jointK(crscor.star, K1, K2)
  KstarKt <- chol(K.big)
  KstarKt.inv <- solve(KstarKt)
  Kstar.inv <- KstarKt.inv %*% t(KstarKt.inv)
  detKstarInv_apx <- -2*sum(log(diag(KstarKt)))
  
  detKInv_apx <- -2*sum(log(diag(Kchol)))
  
  log_rhotStar <-  0.5*detKstarInv_apx - 
    0.5/tau * c( t(c(om[,1], om[,2]) - c(mu1, mu2)) %*% Kstar.inv %*% (c(om[,1], om[,2]) - c(mu1, mu2))) +
    (cor.b - 1)* log(1 - crscor.star) + nu.star - exp(nu.star)
  
  log_rhot <- 0.5*detKInv_apx -
    0.5/tau * c( t(c(om[,1], om[,2]) - c(mu1, mu2)) %*% K.inv %*% (c(om[,1], om[,2]) - c(mu1, mu2))) +
    (cor.b - 1)*log(1 - crscor) + nu - exp(nu)
  
  # Acceptance crieterion
  logR <- log_rhotStar - log_rhot
  xcptProb <- exp(logR)
  
  if(runif(1) <= xcptProb) {
    crscor <- crscor.star 
    Kchol <- KstarKt
    Kchol.inv <- KstarKt.inv
    K.inv <- Kstar.inv
    acceptCountRho <- acceptCountRho + 1
  }
  
  
  # prediction
  k1 <- covar.sep(U.test, U.test, d=-1/(4*log(phi1)), g=nug)
  k2 <- covar.sep(U.test, U.test, d=-1/(4*log(phi2)), g=nug)
  
  k1x <- covar.sep(U.test, U, d=-1/(4*log(phi1)), g=0)
  k2x <- covar.sep(U.test, U, d=-1/(4*log(phi2)), g=0)
  k.joint <- compute.jointK(crscor, k1x, k2x)
  theta.new.mean <- 
    c(rep(beta1, nrow(U.test)), rep(beta2, nrow(U.test))) + 
    k.joint %*% K.inv %*% (logit(c(Theta[,1], Theta[,2])) - c(mu1, mu2))
  theta.new.var <- tau * (compute.jointK(crscor, k1, k2) - k.joint %*% K.inv %*% t(k.joint))
  theta.new.var <- diag(pmax(1e-8, diag(theta.new.var)))
  theta.test <- t(chol(theta.new.var)) %*% rnorm(2*n.test) + theta.new.mean
  test.Theta <- cbind(logit.inv(theta.test[1:nrow(U.test)]), logit.inv(theta.test[(nrow(U.test)+1):(2*nrow(U.test))]))
  
  kx.r <- covar.sep(matrix(x.test/x.max,ncol=1), matrix(x/x.max,ncol=1), d=-1/(4*log(phi.r)), g=0)
  ratio.new.mean <- beta.r + kx.r %*% K.r.inv %*% (logit(ratio) - beta.r)
  ratio.new.var <-  tau.r * (covar.sep(matrix(x.test/x.max,ncol=1), matrix(x.test/x.max,ncol=1), d=-1/(4*log(phi.r)), g=nug) - kx.r %*% K.r.inv %*% t(kx.r))
  ratio.test <- t(chol(ratio.new.var)) %*% rnorm(n.test) + ratio.new.mean
  test.ratio <- logit.inv(ratio.test)
  
  test.Y <- sim.f(c(x,x.test), c(Theta[,1], test.Theta[,1]),
                  c(Theta[,2], test.Theta[,2]))[(n+1):(n+n.test)]
  test.Y <- test.ratio * test.Y
  
  if(i <= nburnin){
    if (i%%100 == 0){
      if (acceptCount/(100*5) > 0.25){ 
        cTheta <- 5.75 * cTheta
      }else if(acceptCount/(100*5) < 0.2){ 
        cTheta <- 0.25 * cTheta
      }
      
      if (acceptCount.ratio/(100*5) > 0.25){ 
        cRatio <- 5.75 * cRatio
      }else if(acceptCount.ratio/(100*5) < 0.2){ 
        cRatio <- 0.25 * cRatio
      }
      
      cPhi1.tmp <- cPhi1
      cPhi1.tmp[acceptCountPhi1/100 > 0.5] <- 1.75 * cPhi1.tmp[acceptCountPhi1/100 > 0.5]
      cPhi1.tmp[acceptCountPhi1/100 < 0.4] <- 0.25 * cPhi1.tmp[acceptCountPhi1/100 < 0.4]
      cPhi1 <- cPhi1.tmp
      
      cPhi2.tmp <- cPhi2
      cPhi2.tmp[acceptCountPhi2/100 > 0.5] <- 1.75 * cPhi2.tmp[acceptCountPhi2/100 > 0.5]
      cPhi2.tmp[acceptCountPhi2/100 < 0.4] <- 0.25 * cPhi2.tmp[acceptCountPhi2/100 < 0.4]
      cPhi2 <- cPhi2.tmp
      
      if (acceptCountPhi.r/100 > 0.5){
        cPhi.r <- 1.75 * cPhi.r
      }else if(acceptCountPhi.r/100 < 0.4){
        cPhi.r <- 0.25 * cPhi.r
      }
      
      if (acceptCountRho/100 > 0.5){
        cRho <- 1.75 * cRho
      }else if(acceptCountRho/100 < 0.4){
        cRho <- 0.25 * cRho
      }
      
      acceptCount <- acceptCount.ratio <- 0
      acceptCountPhi1 <- acceptCountPhi2 <- rep(0, ncol(U))
      acceptCountPhi.r <- 0
      acceptCountRho <- 0
    }
  }else{
    if(i %% 2 == 0){
      Theta.sample[,,(i-nburnin)/2] <- Theta
      ratio.sample[(i-nburnin)/2,] <- ratio
      testTheta.sample[,,(i-nburnin)/2] <- test.Theta
      testRatio.sample[(i-nburnin)/2,] <- test.ratio
      testY.sample[(i-nburnin)/2,] <- test.Y
      tau.sample[(i-nburnin)/2] <- tau
      phi1.sample[(i-nburnin)/2,] <- phi1
      phi2.sample[(i-nburnin)/2,] <- phi2
      phi.r.sample[(i-nburnin)/2] <- phi.r
      beta1.sample[(i-nburnin)/2] <- beta1
      beta2.sample[(i-nburnin)/2] <- beta2
      beta.r.sample[(i-nburnin)/2] <- beta.r
      tau.r.sample[(i-nburnin)/2] <- tau.r
      crscor.sample[(i-nburnin)/2] <- crscor
    }
  }
}

testY.sample <- matrix(rpois(length(testY.sample), testY.sample), ncol = ncol(testY.sample))
