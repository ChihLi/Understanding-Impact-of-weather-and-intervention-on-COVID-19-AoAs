# predicting R0 on the test data (unew)
emulator.fun <- function(unew){
  if(class(unew)[1] == "data.frame")  unew <- as.matrix(unew, nrow=1)
  if(is.null(nrow(unew))) unew <- matrix(unew, nrow=1, ncol=ncol(U))
  
  out <- array(0, c(nrow(unew), n.emulator.MC.sample, 3))
  for(j in 1:n.emulator.MC.sample){
    i <- nsamples/2/n.emulator.MC.sample*(j-1)+1
    
    K1 <- covar.sep(U, d=-1/(4*log(phi1.sample[i,])), g=nug)
    K2 <- covar.sep(U, d=-1/(4*log(phi2.sample[i,])), g=nug)
    K.big <- compute.jointK(crscor.sample[i], K1, K2)
    Kchol <- chol(K.big)
    Kchol.inv <- solve(Kchol)
    K.inv <- Kchol.inv %*% t(Kchol.inv)
    
    k1 <- covar.sep(unew, unew, d=-1/(4*log(phi1.sample[i,])), g=nug)
    k2 <- covar.sep(unew, unew, d=-1/(4*log(phi2.sample[i,])), g=nug)
    
    k1x <- covar.sep(unew, U, d=-1/(4*log(phi1.sample[i,])), g=0)
    k2x <- covar.sep(unew, U, d=-1/(4*log(phi2.sample[i,])), g=0)
    k.joint <- compute.jointK(crscor.sample[i], k1x, k2x)
    theta.new.mean <- 
      c(rep(beta1.sample[i], nrow(unew)), rep(beta2.sample[i], nrow(unew))) + 
      k.joint %*% K.inv %*% (logit(c(Theta.sample[,1,i], Theta.sample[,2,i])) - 
                               c(rep(beta1.sample[i], nrow(U)), rep(beta2.sample[i], nrow(U))))
    kKkt <- rep(0,nrow(k.joint))
    for(jj in 1:nrow(k.joint)) kKkt[jj] <- drop(k.joint[jj,] %*% K.inv %*% k.joint[jj,])
    theta.new.var <- tau.sample[i] * (diag(compute.jointK(crscor.sample[i], k1, k2)) - kKkt)
    theta.new.var <- pmax(1e-8, theta.new.var)
    theta.test <- rnorm(2*nrow(unew), mean = theta.new.mean, sd = sqrt(theta.new.var))
    
    # predictions on beta, gamma, and R0
    theta1.pred <- logit.inv(theta.test[1:nrow(unew)])
    theta2.pred <- logit.inv(theta.test[(nrow(unew)+1):(2*nrow(unew))])      
    theta12.pred <- theta1.pred/theta2.pred
    
    out[,j,1] <- theta12.pred
    out[,j,2] <- theta1.pred
    out[,j,3] <- theta2.pred
  }
  return(out)
}

# compute sensitivity index using MC samples
sensitity.fun <- function(i,j=NULL){
  
  n.sa.MC.sample <- 1000
  
  if(is.null(j)){ # sensitivity index for main effects
    x <- sample(x = U[,i], n.sa.MC.sample, replace = TRUE)
    y <- sample(x = U[,6], n.sa.MC.sample, replace = TRUE)
    
    tmp1 <- matrix(0, ncol = ncol(U), nrow = n.sa.MC.sample)
    for(ii in (1:6)[-i])  tmp1[,ii] <- sample(x = U[,ii], n.sa.MC.sample, replace = TRUE)
    tmp1[,i] <- x
    
    tmp2 <- matrix(0, ncol = ncol(U), nrow = n.sa.MC.sample)
    for(ii in (1:6)[-i])  tmp2[,ii] <- sample(x = U[,ii], n.sa.MC.sample, replace = TRUE)
    tmp2[,i] <- x
  }else{ # sensitivity index for interaction effects
    x <- sample(x = U[,i], n.sa.MC.sample, replace = TRUE)
    y <- sample(x = U[,j], n.sa.MC.sample, replace = TRUE)
    
    tmp1 <- matrix(0, ncol = ncol(U), nrow = n.sa.MC.sample)
    for(ii in (1:6)[-c(i,j)])  tmp1[,ii] <- sample(x = U[,ii], n.sa.MC.sample, replace = TRUE)
    tmp1[,j] <- y
    tmp1[,i] <- x
    
    tmp2 <- matrix(0, ncol = ncol(U), nrow = n.sa.MC.sample)
    for(ii in (1:5)[-i])  tmp2[,ii] <- sample(x = U[,ii], n.sa.MC.sample, replace = TRUE)
    tmp2[,j] <- y
    tmp2[,i] <- x
  }
  Z1 <- emulator.fun(tmp1)[,,1]
  Z2 <- emulator.fun(tmp2)[,,1]
  
  return((apply(Z1 * Z2, 2, mean) - apply(Z1, 2, mean) * apply(Z2, 2, mean))/(apply(Z1^2, 2, mean) - (apply(Z1, 2, mean))^2))
}


# main effect of x on variable i (where i=1,...,5)
main.fun <- function(x, i){

  tmp <- matrix(0, nrow = n.effect.MC.sample * length(x), ncol = ncol(U))
  for(ii in 1:6)  tmp[,ii] <- sample(x = U[,ii], n.effect.MC.sample*length(x), replace = TRUE)
  tmp[,i] <- rep(x, each = n.effect.MC.sample)
  
  emulator.out <- emulator.fun(tmp)
  emulator.out <- array(emulator.out, c(n.effect.MC.sample, length(x), n.emulator.MC.sample, 3))
  
  return(apply(emulator.out, c(2,3,4), median))
}

# interaction effect of x on variables i and 6 (where i=1,...,5)
interaction.fun <- function(x, i){
  
  y <- unique(U[,6])
  tmp <- matrix(0, ncol = ncol(U), nrow = n.effect.MC.sample * length(x) * length(y))
  for(ii in (1:5)[-i])  tmp[,ii] <- sample(x = U[,ii], n.effect.MC.sample*length(x)*length(y), replace = TRUE)
  tmp[,c(i,6)] <- matrix(rep(as.matrix(expand.grid(x,y)), each = n.effect.MC.sample), ncol = 2)
  
  emulator.out <- emulator.fun(tmp)
  emulator.out <- array(emulator.out, c(n.effect.MC.sample, length(x)*length(y), n.emulator.MC.sample, 3))
  
  return(apply(emulator.out, c(2,3,4), median)[,,1])
}