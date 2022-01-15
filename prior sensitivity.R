set.seed(1)
prior.rand <- matrix(runif(n.prior*n.para), ncol=n.para)

for(kkk in 1:n.prior){
  
  # prior
  prior <- function(U){
    
    d <- ncol(U)
    
    beta1.mean <- prior.rand[kkk,1]
    beta2.mean <- prior.rand[kkk,2]
    beta.r.mean <- prior.rand[kkk,3]
    
    beta1.var <- prior.rand[kkk,4]
    beta2.var <- prior.rand[kkk,5]
    beta.r.var <- prior.rand[kkk,6]
    
    tau.shape <- prior.rand[kkk,7]
    tau.rate <- prior.rand[kkk,8]
    
    tau.r.shape <- prior.rand[kkk,9]
    tau.r.rate <- prior.rand[kkk,10]
    
    phi1.b <- phi2.b <- prior.rand[kkk,11]
    phi.r.b <- prior.rand[kkk,12]
    cor.b <- prior.rand[kkk,13]
    
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
  
  set.seed(random.seed[jjj])
  source("realdata_load.R")
  source("posterior_sampling.R")
  save(testTheta.sample, file = paste0(metro.city, "_posterior_", kkk,".Rdata"))
}
