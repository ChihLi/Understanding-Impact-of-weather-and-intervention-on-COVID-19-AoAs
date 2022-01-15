# load necessary packages
library(plgp)
library(ggplot2)
library(ggpubr)
library(gbm)
library(gridExtra)
library(PerformanceAnalytics)

t1 <- Sys.time()

#####     reproducing the results in Section 4      ######
## "simulation_study.R" produces three pdf files:       ## 
## - "simulation_setup.pdf" (Fig. 1)                    ##
## - "simulation_correlation.pdf" (Fig. 2)              ##
## - "simulation_prediction.pdf" (Fig. 3)               ##
##########################################################
source("simulation_study.R")

#####      reproducing the results in Section 5      ##### 
## the following command lines produce eight pdf files: ## 
## - "corplot.pdf" (Fig. 4)                             ##
## - "mean_effect.pdf" (Fig. 5)                         ##
## - "main_effect.pdf" (Fig. 6)                         ##
## - "main_effect_sa.pdf" (Fig. 7)                      ##
## - "interaction_plot.pdf" (Fig. 8)                    ##
## - "newcase_prediction.pdf" (Fig. 9)                  ##
## - "realdata_fraction.pdf" (Fig. 10)                  ##
## - "deviation.pdf" (Fig. 11)                          ##
## NOTE: some temporary .RData files will be            ##
##       produced along these pdf files.                ##
##       These temporary files will be removed          ##
##       after running the last command line            ##                  
##########################################################
source("cor.plot.R")
source("logit.R")
source("init.R")
source("prior.R")
source("SIR.R")
source("sensitivity_funcs.R")

### read US county level data
# specify the cities and states
metro.city.all <- c("New York", "Los Angeles", "Baltimore", "Houston", "San Francisco", "Chicago", "Saint Louis", "Atlanta")
metro.state.all <- c("New York", "California", "Maryland", "Texas", "California", "Illinois", "Missouri", "Georgia")
n.city <- length(metro.city.all)
n.test <- 14 #from Nov 12nd to 25th, 2020
all.data <- read.csv("data/us-counties.csv")
county.name <- read.delim("data/county.txt", header = FALSE, stringsAsFactors = FALSE)
county.name <- county.name[,c(3,4,6)]
county.name <- unique(county.name)
colnames(county.name) <- c("city", "state", "county")

# setup for MCMC
nsamples <- 2000
nburnin <- 2000
n.emulator.MC.sample <- 100
n.sa.MC.sample <- 1000
n.effect.MC.sample <- 100
n.main.rand <- 30 # number of uniform quantiles of x for main effects
n.int.rand <- 20 # number of uniform quantiles of x for interaction effects
random.seed <- c(10, 20, 300, 40, 5, 6000, 7, 80) # could try different chains

# number of the prior specifications for sensitivity analysis for priors
n.prior <- 30 
n.para <- 13

# this can take a long time
# if parallel computing or submitting jobs to cluster is available, 
# I would suggest to do so (for each jjj).
for(jjj in 1:n.city){ 
  
  metro.city <- metro.city.all[jjj]
  metro.state <- metro.state.all[jjj]
  cat("Running", metro.city, "...\n")
  
  set.seed(random.seed[jjj])
  
  ### load data ###
  cat("=== loading data\n")
  source("realdata_load.R")
  
  ### draw posteriors ###
  cat("=== draw posteriors\n")
  source("posterior_sampling.R")
  # save files
  save(x, y, x.test, y.test, date, y.fit, y.pred,
       testTheta.sample, ratio.sample, file = paste0(metro.city, "_posterior.Rdata"))
  
  ### perform sensitivity analysis ###
  cat("=== sensitivity analysis\n")
  source("sensitivity analysis.R") 
  
  ### sensitivity analysis for the priors ###
  cat("=== prior sensitivity\n")
  source("prior sensitivity.R")
}

### Reproduce the figures

#####  Fig. 4: scatter plot ##### 
rawdata <- read.csv("data/climate_city.csv")[,-c(1,2,4)]
colnames(rawdata)[6] <- "Intervention" 
pdf("corplot.pdf", width = 6, height = 6)
cor.plot(rawdata, histogram=TRUE)
dev.off()

##### Fig. 5: overall mean effect ##### 
mean.df <- data.frame("city" = rep(metro.city.all, each = n.emulator.MC.sample))
mean.df <- cbind(mean.df, "R0" = NA)
mean.df$city <- factor(mean.df$city , levels=metro.city.all)
for(i in 1:n.city){
  load(paste0(metro.city.all[i],"_mean.Rdata"))
  mean.df[((i-1)*n.emulator.MC.sample+1):(i*n.emulator.MC.sample),2] <- mean.effect[,1]
}
pdf("mean_effect.pdf", width = 5, height = 3)
p <- ggplot(mean.df, aes(x=city, y=R0)) + geom_boxplot()
p + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) + ylab("Mean effect of R0")
dev.off()

##### Fig. 6: main effect sensitivity index ##### 
main.sa.df <- data.frame(matrix(0, nrow=length(metro.city.all) * n.emulator.MC.sample, ncol = 6))
U.name <- c("Temperature", "Humidity", "WindSpeed", "Pressure", "Precipitation", "Intervention")
colnames(main.sa.df) <- c(U.name)
for(i in 1:length(metro.city.all)){
  load(paste0(metro.city.all[i],"_main_sensitivity.Rdata"))
  main.sa.df[((i-1)*n.emulator.MC.sample+1):(i*n.emulator.MC.sample),] <- sensitity.mx
}
main.sa.df <- stack(main.sa.df)
main.sa.df$city <- rep(rep(metro.city.all, each = n.emulator.MC.sample), 6)
main.sa.df$city <- factor(main.sa.df$city , levels=metro.city.all)

pdf("main_effect_sa.pdf", width = 8, height = 5)
m1 <- ggplot(main.sa.df[is.element(main.sa.df$city, metro.city.all[1:4]),], aes(x=ind, y=values) ) + 
  geom_boxplot(alpha = 0.5, show.legend = FALSE) + facet_grid(.~city) + ylab("Sensitivity index") + 
  theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 8)) + ylim(c(-0.1, 0.6))

m2 <- ggplot(main.sa.df[is.element(main.sa.df$city, metro.city.all[5:8]),], aes(x=ind, y=values) ) + 
  geom_boxplot(alpha = 0.5, show.legend = FALSE) + facet_grid(.~city) + ylab("Sensitivity index") + 
  theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 8)) + ylim(c(-0.1, 0.6))
ggarrange(m1,m2, ncol=1, nrow=2)
dev.off()


##### Fig. 7: plot main effects ##### 
main.df <- data.frame("city" = factor(rep(metro.city.all, each = n.main.rand)), matrix(NA,ncol=12,nrow=n.main.rand * length(metro.city.all)))
colnames(main.df)[-1] <- c(U.name, paste0(U.name,".R0"))
main.df$city <- factor(main.df$city , levels=metro.city.all)

for(i in 1:length(metro.city.all)){
  load( paste0(metro.city.all[i],"_main_effect.Rdata"))
  main.df[ ((i-1)*n.main.rand+1):(n.main.rand*i),U.name] <- U.rand
  for(ii in 1:6) {
    main.df[((i-1)*n.main.rand+1):((i-1)*n.main.rand+nrow(R0.main.ls[[ii]])),paste0(U.name[ii],".R0")] <-  apply(R0.main.ls[[ii]],1,median)
  }
}

pdf("main_effect.pdf", width = 6, height = 6.5)
theme_set(theme_bw())
r1 <- ggplot(main.df, aes(x=Temperature, y=Temperature.R0, group=city)) +
  geom_line(aes(linetype=city, size=city, colour=city)) + xlim(c(40,80)) + coord_cartesian(ylim=c(-0.2,0.3)) + 
  scale_colour_manual(values=paste0("gray",(0:7)*10))+
  scale_size_manual(values=seq(0.3,1,length.out = 8)) +
  theme(legend.position = "none", axis.title.y=element_text(size=8)) + xlab("Temperature (°F)")+ ylab("R0 - mean effct of R0") 
r2 <- ggplot(main.df, aes(x=Humidity, y=Humidity.R0, group=city)) +
  geom_line(aes(linetype=city, size=city, colour=city)) + xlim(c(50,85)) + coord_cartesian(ylim=c(-0.2,0.3)) + 
  scale_colour_manual(values=paste0("gray",(0:7)*10))+
  scale_size_manual(values=seq(0.3,1,length.out = 8)) +
  theme(legend.position = "none", axis.title.y=element_text(size=8)) + xlab("Humidity (%)")+ ylab("R0 - mean effct of R0") 
r3 <- ggplot(main.df, aes(x=WindSpeed, y=WindSpeed.R0, group=city)) +
  geom_line(aes(linetype=city, size=city, colour=city)) + xlim(c(2.5,15)) + coord_cartesian(ylim=c(-0.2,0.3)) + 
  scale_colour_manual(values=paste0("gray",(0:7)*10))+
  scale_size_manual(values=seq(0.3,1,length.out = 8)) +
  theme(legend.position = "none", axis.title.y=element_text(size=8)) + xlab("Wind Speed (mph)")+ ylab("R0 - mean effct of R0") 
r4 <- ggplot(main.df, aes(x=Pressure, y=Pressure.R0, group=city)) +
  geom_line(aes(linetype=city, size=city, colour=city)) + coord_cartesian(ylim=c(-0.2,0.3)) + 
  scale_colour_manual(values=paste0("gray",(0:7)*10))+
  scale_size_manual(values=seq(0.3,1,length.out = 8)) +
  theme(legend.position = "none", axis.title.y=element_text(size=8)) + xlab("Pressure (Hg)")+ ylab("R0 - mean effct of R0") 
r5 <- ggplot(main.df, aes(x=Precipitation, y=Precipitation.R0, group=city)) +
  geom_line(aes(linetype=city, size=city, colour=city))+ xlim(c(0,1)) + coord_cartesian(ylim=c(-0.2,0.3))   + 
  scale_colour_manual(values=paste0("gray",(0:7)*10))+
  scale_size_manual(values=seq(0.3,1,length.out = 8)) +
  theme(legend.position = "none", axis.title.y=element_text(size=8)) + xlab("Precipitation (inches)")+ ylab("R0 - mean effct of R0") 
r6 <- ggplot(main.df, aes(x=Intervention, y=Intervention.R0, group=city)) +
  geom_line(aes(linetype=city, size=city, colour=city)) +
  scale_colour_manual(values=paste0("gray",(0:7)*10))+
  scale_size_manual(values=seq(0.3,1,length.out = 8)) + coord_cartesian(ylim=c(-0.2,0.3))+ 
  theme(legend.position = "none", axis.title.y=element_text(size=8)) + xlab("Intervention level")+ ylab("R0 - mean effct of R0") 

ggarrange(r1, r2, r3, r4, r5, r6, ncol=2, nrow=3, common.legend = TRUE, legend="top")
dev.off()


##### Fig. 8: plot interaction effects ##### 
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

##### Fig. 9: predictions on the test data ##### 
pdf("newcase_prediction.pdf", width = 8, height = 8)
theme_set(theme_bw())
gg2.ls <- vector("list", 8)
for(jjj in 1:8){
  metro.city <- metro.city.all[jjj]
  load(paste0(metro.city, "_posterior.Rdata"))
  
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


##### Fig. 10: fractions of cases reported ##### 
pdf("realdata_fraction.pdf", width = 6, height = 6)
layout(matrix(1:n.city,ncol=2,byrow=TRUE))
for(jjj in 1:n.city){
  metro.city <- metro.city.all[jjj]
  load(paste0(metro.city, "_posterior.Rdata"))
  fraction.df <- data.frame(date=date[1:length(x)],
                            ratio=apply(ratio.sample,2,mean),
                            y=y[1:length(x)])
  sample.df <- data.frame(date=date[1:length(x)],
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

##### Fig. 11: sensitivity analysis for priors ##### 
dev.R0.out <- matrix(0, nrow=n.test*n.prior, ncol=n.city)
colnames(dev.R0.out) <- metro.city.all
for (metro.city in metro.city.all) {
  load(paste0(metro.city, "_posterior.Rdata"))
  posterior.R0.ori <- apply(testTheta.sample,c(1,2),median)
  posterior.R0.ori <- posterior.R0.ori[,1]/posterior.R0.ori[,2]
  
  posterior.R0.new <- matrix(0,nrow=length(posterior.R0.ori),ncol=n.prior)
  for(kkk in 1:n.prior){
    load(paste0(metro.city,"_posterior_", kkk,".Rdata"))
    Theta.tmp <- apply(testTheta.sample,c(1,2),median)
    posterior.R0.new[,kkk] <- Theta.tmp[,1]/Theta.tmp[,2]
  }
  
  dev.R0.out[,metro.city] <- c(apply(posterior.R0.new, 2, function(x) (x-posterior.R0.ori)/posterior.R0.ori))
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

### remove all .RData files ###
all.files <- list.files()
rdata.files <- grep(".Rdata", all.files, value = TRUE)
file.remove(rdata.files)

t2 <- Sys.time()
cat("it takes", difftime(t2,t1,units="days"), "days")
