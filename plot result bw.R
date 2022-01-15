##### interaction effect ##### 
jjj <- 1 # run New York
load(paste0(metro.city.all[jjj],"_interaction_effect.Rdata"))
interac.df <- data.frame(matrix(0, nrow=100, ncol = 4))
colnames(interac.df) <- c("X", "level", "U", "city")

interac.df <- data.frame(X = c(interaction.effect.mx), 
                         level = rep(round((U6.unique * scale.U[6] + center.U[6])^2), each = nrow(interaction.effect.mx)),
                         U = rep(U.rand, length(U6.unique)))
interac.df$level <- as.factor(interac.df$level)

load(paste0(metro.city.all[jjj],"_interaction_sensitivity.Rdata"))
colnames(sensitity.mx) <- U.name[-6]
interact.sa.df <- stack(data.frame(sensitity.mx))

pdf("interaction_plot.pdf", width = 8, height = 3)
g1 <- ggplot(interact.sa.df, aes(x=ind, y=values) ) + 
  geom_boxplot(alpha = 0.5, show.legend = FALSE) + ylab("Sensitivity index") + 
  theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 8))

g2 <- ggplot(interac.df, aes(x=U, y=X, linetype=level, color=level) ) + 
  geom_line(aes(size=level)) + scale_size_manual(values=seq(0.2,1,0.2)) + 
  scale_colour_manual(values=paste0("gray",(0:4)*10))+
  ylab("Sensitivity index") + labs(x = "Temperature", y = "R0") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(5.5, 5.5, 30.5, 5.5, "pt"))

grid.arrange(g1, g2, ncol = 2, widths = c(10,12))
dev.off()


##### predictions on the test data ##### 
pdf("newcase_prediction.pdf", width = 8, height = 8)
theme_set(theme_bw())
gg2.ls <- vector("list", 8)
for(jjj in 1:8){
  source("realdata_load.R")
  load(paste0(metro.city, "_posterior.RData"))
  date <- subdata[-1,"date"]
  y.fit <- matrix(0, nrow=length(x), ncol=1000)
  for(i in 1:1000)  y.fit[,i] <- ratio.sample[i,] * sim.f(x, Theta.sample[,1,i], Theta.sample[,2,i], city)
  y.test[y.test<0] <- 0 
  y.pred <- matrix(rpois(length(testY.sample), testY.sample), ncol = ncol(testY.sample))
  
  data.df <- data.frame(date=c(date[1:length(y)], date[length(y)+1], date[(length(y)+1):length(date)]), 
                        y=c(y,y.test[1],y.test), data=c(rep("train",length(y)+1), rep("test",length(y.test))),
                        pred = c(apply(y.fit, 1, mean), mean(y.pred[,1]), apply(y.pred, 2, mean)),
                        lb = c(rep(NA,length(y)+1),apply(y.pred, 2, min)), 
                        up = c(rep(NA,length(y)+1),apply(y.pred, 2, max)))
  
  data.df <- data.df[data.df$date >= "2020-07-01", ]
  gg2 <- ggplot(data.df, aes(x=date, y= pred,linetype=data),colour="black")  + 
    geom_point(aes(x=date, y= y), shape = 21, colour = "gray45", fill = "white", size = 1.5) + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous(labels = scales::comma) + 
    scale_linetype_manual(values=c("solid", "dashed"))+
    labs( y = "daily cases", title = metro.city) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))
  
  gg2.ls[[jjj]] <- gg2 + geom_line(aes(size=data)) + scale_size_manual(values=c(1,0.5))+ geom_smooth(aes(ymin = lb, ymax = up), stat = "identity",colour="black",size=0.6) +  theme(legend.position = "none")
}
grid.arrange(gg2.ls[[1]], gg2.ls[[2]], gg2.ls[[3]], gg2.ls[[4]], 
             gg2.ls[[5]], gg2.ls[[6]], gg2.ls[[7]], gg2.ls[[8]], ncol = 2)
dev.off()


##### fractions of cases reported ##### 
pdf("realdata_fraction.pdf", width = 6, height = 6)
layout(matrix(1:8,ncol=2,byrow=TRUE))
for(jjj in 1:8){
  source("realdata_load.R")
  load(paste0(metro.city, "_posterior.RData"))
  fraction.df <- data.frame(date=subdata[-1,"date"][1:length(x)],
                            ratio=apply(ratio.sample,2,mean),
                            y=y[1:length(x)])
  sample.df <- data.frame(date=subdata[-1,"date"][1:length(x)],
                          upper=apply(ratio.sample,2,quantile,0.975),
                          lower=apply(ratio.sample,2,quantile,0.025))
  par(mgp=c(3,1,0))
  par(mar=c(3, 4, 1.5, 4) + 0.1)
  plot(fraction.df[,c("date","y")],xlab="", ylab="daily cases", las=2, col="gray", cex=0.6, cex.axis=0.9)
  mtext(side = 3, line = 0.5, cex=0.8, metro.city)
  par(new=TRUE)
  plot(fraction.df[,c("date","ratio")],type="n",axes=FALSE,ylab="",xlab="",ylim=c(0,1))
  lines(sample.df[,c("date","upper")],lty=2)
  lines(sample.df[,c("date","lower")],lty=2)
  lines(fraction.df[,c("date","ratio")],lwd=2.5)
  axis(4, ylim=c(0,1),las=1, cex.axis=0.9)
  mtext(side = 4, line = 2.5, cex=0.7, "fraction")
}
dev.off()

##### sensitivity analysis for priors ##### 
dev.R0.out <- matrix(0, nrow=14*n.prior, ncol=8)
colnames(dev.R0.out) <- metro.city.all
for (metro.city in metro.city.all) {
  load(paste0(metro.city, "_posterior.RData"))
  posterior.R0.ori <- apply(testTheta.sample,c(1,2),median)
  posterior.R0.ori <- posterior.R0.ori[,1]/posterior.R0.ori[,2]

  posterior.R0.new <- matrix(0,nrow=length(posterior.R0.ori),ncol=n.prior)
  for(kkk in 1:n.prior){
    load(paste0(metro.city,"_posterior_", kkk,".RData"))
    Theta.tmp <- apply(testTheta.sample,c(1,2),median)
    posterior.R0.new[,kkk] <- Theta.tmp[,1]/Theta.tmp[,2]
  }
  
  dev.R0.out[,metro.city] <- apply(posterior.R0.new, 2, function(x) (x-posterior.R0.ori)/posterior.R0.ori)
}

dev.R0.boxplot <- stack(data.frame(dev.R0.out))
colnames(dev.R0.boxplot) <- c("deviation", "city")
dev.R0.boxplot[,1] <- dev.R0.boxplot[,1] * 100
dev.R0.boxplot[,2] <- sub("\\.", " ", dev.R0.boxplot[,2])
dev.R0.boxplot$city <- factor(dev.R0.boxplot$city , levels=metro.city.all)

pdf("deviation.pdf", width = 5, height = 3)
p <- ggplot(dev.R0.boxplot, aes(x=city, y=deviation)) + geom_boxplot()
p + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) + 
  ylab("Percentage deviation (%)")+ coord_cartesian(ylim=c(-20,80)) 
dev.off()
