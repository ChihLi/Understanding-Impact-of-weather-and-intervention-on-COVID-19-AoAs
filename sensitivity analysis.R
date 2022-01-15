# settings
set.seed(random.seed[jjj])

# sensitivity index for mean effect
X.rand <- matrix(0, nrow = n.sa.MC.sample, ncol = ncol(U))
for(ii in 1:6) X.rand[,ii] <- sample(x = U[,ii], n.sa.MC.sample, replace = TRUE)
mean.effect <- apply(emulator.fun(X.rand), c(2,3), median, na.rm=TRUE)
save(mean.effect, file = paste0(metro.city,"_mean.Rdata"))
  
# sensitivity index for main effects
sensitity.mx <- matrix(0, nrow = n.emulator.MC.sample, ncol=6)
for(i in 1:6) sensitity.mx[,i] <- sensitity.fun(i)
save(sensitity.mx, file = paste0(metro.city,"_main_sensitivity.Rdata"))

# main effects
U.rand <- matrix(NA, nrow = n.main.rand, ncol = ncol(U))
R0.main.ls <- vector("list", ncol(U)) 
for(i in 1:ncol(U)){
  if(i == 5){
    x.rand <- sort(unique(U[,i]))
    if(length(x.rand) > n.main.rand){
      x.rand <- x.rand[seq(1, length(x.rand), length.out = n.main.rand)]
    }
    U.rand[1:length(x.rand),i] <- (x.rand * scale.U[i] + center.U[i])^2  - 1
  }else if(i == 6){
    x.rand <- sort(unique(U[,i]))
    U.rand[1:length(x.rand),i] <- (x.rand * scale.U[i] + center.U[i])^2
  }else{
    # x.rand: uniform quantiles
    x.rand <- seq(quantile(U[,i],0.05),quantile(U[,i],0.95), length.out = n.main.rand)
    U.rand[1:length(x.rand),i] <- x.rand * scale.U[i] + center.U[i]
  }
  main.effect.ls <- main.fun(x.rand, i)
  R0.main.ls[[i]] <- t(t(main.effect.ls[,,1]) - mean.effect[,1])
}
save(U.rand, R0.main.ls, file = paste0(metro.city,"_main_effect.Rdata"))

# sensitivity for interaction effects
sensitity.mx <- matrix(0, nrow = n.emulator.MC.sample, ncol=5)
for(i in 1:5) sensitity.mx[,i] <- sensitity.fun(i,6)
save(sensitity.mx, file = paste0(metro.city,"_interaction_sensitivity.Rdata"))

# interaction effects
R0.int.ls <- vector("list", ncol(U))
U6.unique <- unique(U[,6])
plot.index <- which.max(apply(sensitity.mx, 2, mean)) # only plot the strongest effect
x.rand <- seq(quantile(U[,plot.index],0.05),quantile(U[,plot.index],0.95), length.out = n.int.rand)
if(plot.index == 5) U.rand <- (x.rand * scale.U[plot.index] + center.U[plot.index])^2-1 else U.rand <- x.rand * scale.U[plot.index] + center.U[plot.index]
interaction.effect <- interaction.fun(x.rand, plot.index)
interaction.effect <- apply(interaction.effect, 1, mean)
interaction.effect.mx <- matrix(interaction.effect, nrow= n.int.rand)
save(U.rand, U6.unique, scale.U, center.U, interaction.effect.mx, file = paste0(metro.city,"_interaction_effect.Rdata"))  

