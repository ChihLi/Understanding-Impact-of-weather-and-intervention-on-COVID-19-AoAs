### specify county
if(any(metro.city == "San Francisco")){
  county.select <- county.name[is.element(county.name[,"city"], c(metro.city, "San Jose"))
                               & is.element(county.name[,"state"], metro.state),]
  county.select$city[county.select$city == "San Jose"] <- "San Francisco"
}else{
  county.select <- county.name[is.element(county.name[,"city"], metro.city) 
                               & is.element(county.name[,"state"], metro.state),]
  county.select$county[county.select$county == "New York"] <- "New York City"
  county.select$county[county.select$county == "Baltimore (city)"] <- "Baltimore city"
}

### find the sub-data for metro.city's county
select.data <- all.data[is.element(all.data[,"county"], county.select$county) & 
                          is.element(all.data[,"state"], county.select$state),]
colnames(select.data)[colnames(select.data) == "county"] <- "city"
for(i in 1:nrow(county.select)){
  select.data$city <- gsub(county.select$county[i], county.select$city[i], select.data$city)
}
select.data$date <- as.Date(select.data$date,"%m/%d/%y")

select.data <- aggregate(select.data[,c("cases", "deaths")],
                         list("date" = select.data$date, "city" = select.data$city, "state" = select.data$state), sum)

### read population
population.all <- read.csv("data/city_population.csv")

# y is daily case and x is date
n.days <- table(select.data[,"city"])[metro.city] - 1
y <- rep(0, sum(n.days) - n.test)
x <- rep(0, sum(n.days) - n.test)
y.test <- rep(0, n.test)
x.test <- rep(0, n.test)

population <- population.all[population.all$city == metro.city, "population"]
subdata <- select.data[select.data$city == metro.city, ]
date <- subdata[-1,"date"]
y.all <- diff(subdata$cases)
if(metro.city == "Atlanta"){ # abnormal data
  y.all[217] <- round(mean(y.all[214:220]))
}

if(metro.city == "Houston"){ # abnormal data
  y.all[200] <- round(mean(y.all[194:206]))
}
y <- y.all[1:(length(y.all)-n.test)]
x <- 1:length(y)

y.test <- y.all[(length(y.all)-n.test+1):length(y.all)]
x.test <- (length(y)+1):(length(y)+n.test)

I0 <- subdata$cases[1]
R0 <- subdata$deaths[1]
S0 <- population - I0 - R0
N0 <- population

first.incidence <- as.character(subdata$date[1])

### training input data: weather and government interventions
U <- read.csv("data/climate_city.csv")
U <- U[is.element(U$city, metro.city),]
U <- as.matrix(U[,c(3,5,6,7,8,9),drop=FALSE])
U <- U[-1,]
U[,5] <- (U[,5]+1)^{1/2}
U[,6] <- sqrt(U[,6])

# scaling the inputs 
U <- scale(U)/sqrt(12)
center.U <-  attr(U,"scaled:center")
scale.U <-  attr(U,"scaled:scale")*sqrt(12)

### test input data: weather and government interventions
U.test <- read.csv("data/climate_city_test.csv")
U.test <- U.test[is.element(U.test$city, metro.city),]
U.test <- as.matrix(U.test[,c(3,5,6,7,8,9),drop=FALSE])
U.test[,5] <- (U.test[,5]+1)^{1/2}
U.test[,6] <- sqrt(U.test[,6])

# scaling the inputs 
U.test <- t((t(U.test)-center.U)/scale.U)
