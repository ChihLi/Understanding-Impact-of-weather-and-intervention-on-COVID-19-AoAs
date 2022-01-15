set.seed(1) # for reproduction
n <- 80
n.train <- 60
n.test <- n - n.train
x.all <- 1:n
U.all <- matrix(seq(0, 1, length.out = n)[sample(n)], ncol =1)

true.mean <- function(U, x){
  exp(-x/50)*(5*(sin(3*U)*exp(-U)+0.2)+ sin(3*U)*(x/10)^2+(x/10))^2
}

y.all <- rpois(n, true.mean(U.all, x.all))
U <-  U.all[1:n.train,,drop=FALSE]
U.test <-  U.all[(n.train+1):n,,drop=FALSE]
x <- x.all[1:n.train]
x.test <- x.all[(n.train+1):n]
y <- y.all[1:n.train]
y.test <- y.all[(n.train+1):n]

sim.f <- function(x, theta1, theta2){
  (5*theta1+theta2*(x/10)^2+(x/10))^2
}

pdf("simulation_setup.pdf", width = 7, height = 2.3)

par(mfrow=c(1,3), mar = c(4, 4, 1, 0) + 0.3, mgp = c(3, 1.5, 0))  
curve(sin(x*3)*exp(-x)+0.2, ylim=c(0,1),xlab="x",ylab=expression(paste(beta(x), " ,",gamma(x))))
curve(sin(x*3),add=TRUE,lty =2)
curve(exp(-x/50),xlab="t",ylab=expression(kappa(t)), xlim=c(0,n))
plot(true.mean(U.all, x.all), xlab = "t",type="l", ylab = "y(t)", lty=2)
points(y.all, cex=1.5)
dev.off()

source("simulation_jointGP.R")

pdf("simulation_correlation.pdf", width = 7, height = 2.3)
par(mfrow=c(1,3), mar = c(4, 4, 1, 0) + 0.3, mgp = c(3, 1.5, 0))  
curve(sin(x*3)*exp(-x)+0.2, ylim=c(0,1),xlab="x",ylab=expression(beta(x)),type="n")
for(j in 1:1000)  lines(sort(U.all), c(Theta.sample[,1,j], testTheta.sample[,1,j])[sort.int(U.all,index.return = TRUE)$ix], col = "gray")
lines(sort(U.all), c(apply(Theta.sample[,1,],1, median), apply(testTheta.sample[,1,],1, median))[sort.int(U.all,index.return = TRUE)$ix], col = 1, lwd = 2)
curve(sin(x*3)*exp(-x)+0.2, ylim=c(0,1),
      xlab="x",ylab=expression(beta(x)), lty = 2, add = TRUE, lwd = 1)

curve(sin(x*3),xlab="x",ylab=expression(gamma(x)),type="n")
for(j in 1:1000)  lines(sort(U.all), c(Theta.sample[,2,j], testTheta.sample[,2,j])[sort.int(U.all,index.return = TRUE)$ix], col = "gray")
lines(sort(U.all), c(apply(Theta.sample[,2,],1, median), apply(testTheta.sample[,2,],1, median))[sort.int(U.all,index.return = TRUE)$ix], col = 1, lwd = 2)
curve(sin(x*3), ylim=c(0,1),
      xlab="x",ylab=expression(gamma(x)), lty = 2, add = TRUE, lwd = 1)

curve(exp(-x/50),xlab="t",ylab=expression(kappa(t)), xlim=c(0,n),type="n")
for(j in 1:1000)  lines(x.all, c(ratio.sample[j,], testRatio.sample[j,]), col = "gray")
lines(x.all, c(apply(ratio.sample,2, median), apply(testRatio.sample,2, median)), col = 1, lwd = 2)
curve(exp(-x/50),
      xlab="x",ylab=expression(kappa(t)), lty = 2, add = TRUE, lwd = 1)
dev.off()

pdf("simulation_prediction.pdf", width = 6, height = 3)
par(mfrow=c(1,1), mar = c(4, 4, 1, 0) + 0.3, mgp = c(3, 1.5, 0))  
ts.plot(true.mean(U.all, x.all),
        xlab = "t", ylab = "y(t)", lty = 2, lwd = 1, ylim=c(0,1400))
for(j in 1:1000){
  lines(x.test, testY.sample[j,], col = "gray")
}
lines(x.test, apply(testY.sample, 2, median), col = 1, lwd = 2)
lines(x.test, true.mean(U.test, x.test), col = 1, lty = 2, lwd = 1)
points(y.all, cex=1.5)
dev.off()
