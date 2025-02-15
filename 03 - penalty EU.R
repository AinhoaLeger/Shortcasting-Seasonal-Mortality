
######################################################
##########     PENALTY ORDER FOR RATES     ##########
######################################################
#####-------------- Directory and packages --------------#####

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
load("fig/data/series.EU.RData")

# devtools::install_github("teunbrand/ggh4x"force = TRUE)
library(tidyverse)
require(ggplot2)
library(splines)
library(zoo)
library(magic)
library(geofacet)
library(ggh4x)
library(colorspace)
library(knitr)
library(kableExtra)

# Load function for forecasting
source("fun/ForFun total.R")

#####-------------- Function for rates --------------#####

period <- sort(unique(death$date))

# Names of the sub-population and accuracy measures
accuracy.names <- c("bic.t1s1","bic.t2s2","bic.t1s2","bic.t2s1",
                    "rmse.t1s1","rmse.t2s2","rmse.t1s2","rmse.t2s1",
                    "mape.t1s1","mape.t2s2","mape.t1s2","mape.t2s1")
country.names= c("Bulgaria","Czechia","Denmark","Germany",
                 "Estonia","Ireland","Greece","Spain","France",
                 "Croatia","Italy","Lithuania","Luxembourg",
                 "Hungary","Netherlands","Austria","Poland",
                 "Portugal","Romania","Slovenia","Finland",
                 "Sweden","Iceland","Norway","Switzerland")

accuracy <- function(data, country.names, country.db) {
  
  # Initialize the dataset with bic, rmse, mape
  bic.t1s1 <- bic.t2s2 <- bic.t1s2 <- bic.t2s1 <-
    rmse.t1s1 <- rmse.t2s2 <- rmse.t1s2 <- rmse.t2s1 <-
    mape.t1s1 <- mape.t2s2 <- mape.t1s2 <- mape.t2s1 <-  c()
  
  for(j in 1:length(country.names)){
    for(i in window){ 
      
      d <- (data %>% filter(country==country.names[j]))$deaths[i:(i+ttot-1)]
      ex <- (data %>% filter(country==country.names[j]))$exposure[i:(i+ttot-1)]
      x <- 1:length(d)
      d.fit <- d[w==1] # subset estimation dataset
      ex.fit <- ex[w==1]
      n.for <- 12 # we produce 12 forecast in a year
      r.obs <- (d/ex)[(tfit+1):ttot] # subset prediction dataset
      
      #####----------------- Grid search for lambdas -----------------######
      
      ##### Grid search for lambda1 and lambda2 (MM) and lambda (SM)
      lambda.seq <- 10^seq(5,7,length.out=30)
      mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
      rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
      mape.vec <- rep(NA, length(lambda.seq)) 
      mape.s.MM <- mape.s.SM <- rep(NA,300)
      
      ##### Rolling window on the estimation set
      for (h in 1:length(lambda.seq)){
        for(g in 1:length(lambda.seq)){
          for(k in c(73,85,97,109)){ # 5 years: 25,37,49
            d.sel <- d.fit[1:(k+11)]
            ex.sel <- ex.fit[1:(k+11)]
            w.sel <- rep(1,k+11)
            w.sel[k:(k+11)] <- 0
            nx.sel <- (floor(k/12)+1)*2
            Output.MM <- PIRLS(d.sel, ex.sel, w.sel, nx.sel, bdeg, 2, 1, 
                               lambda.seq[h], lambda.seq[g])$theta
            Output.SM <- PIRLS2(d.sel, ex.sel, w.sel, nx.sel, bdeg, 2,
                                lambda.seq[g])$theta
            x.s <- 1:length(d.sel)
            X.s <- bspline(x.s, nx.sel, bdeg)
            eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
            eta.s.SM <- c(X.s$BB2 %*% Output.SM)[k:(k+11)]
            r.obs.s <- (d.fit/ex.fit)[k:(k+11)]
            mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
            mape.s.SM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.SM))/r.obs.s)))
          }
          mape.mat[h,g] <- mean(na.omit(mape.s.MM))
          mape.vec[g] <- mean(na.omit(mape.s.SM))
        }}
      mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
      lambda1 = lambda.seq[mape.min.MM[1]]
      lambda2 = lambda.seq[mape.min.MM[2]]
      mape.min.SM <- which(mape.vec == min(mape.vec), arr.ind = TRUE)
      lambda = lambda.seq[mape.min.SM[1]]
      
      #####----------------- Modulation model -----------------######
      
      # Modulation model t1s1
      theta.MM <- PIRLS(d, ex, w, nx, bdeg, 1, 1, lambda1, lambda2)
      X.MM <- bspline(x, nx, bdeg)$BB
      for.t1s1 <- c(X.MM %*% theta.MM$theta)[(tfit+1):ttot]
      bic.t1s1[i] <- theta.MM$Bic
      rmse.t1s1[i] <- (sqrt((1/n.for)*(sum((exp(for.t1s1)-r.obs)^2))))
      mape.t1s1[i] <- ((100/n.for)*sum(abs((r.obs-exp(for.t1s1))/r.obs)))
      
      # Modulation model t2s2
      theta.MM <- PIRLS(d, ex, w, nx, bdeg, 2, 2, lambda1, lambda2)
      X.MM <- bspline(x, nx, bdeg)$BB
      for.t2s2 <- c(X.MM %*% theta.MM$theta)[(tfit+1):ttot]
      bic.t2s2[i] <- theta.MM$Bic
      rmse.t2s2[i] <- (sqrt((1/n.for)*(sum((exp(for.t2s2)-r.obs)^2))))
      mape.t2s2[i] <- ((100/n.for)*sum(abs((r.obs-exp(for.t2s2))/r.obs)))      
      
      # Modulation model t1s2
      theta.MM <- PIRLS(d, ex, w, nx, bdeg, 1, 2, lambda1, lambda2)
      X.MM <- bspline(x, nx, bdeg)$BB
      for.t1s2 <- c(X.MM %*% theta.MM$theta)[(tfit+1):ttot]
      bic.t1s2[i] <- theta.MM$Bic
      rmse.t1s2[i] <- (sqrt((1/n.for)*(sum((exp(for.t1s2)-r.obs)^2))))
      mape.t1s2[i] <- ((100/n.for)*sum(abs((r.obs-exp(for.t1s2))/r.obs)))      
      
      # Modulation model t2s1
      theta.MM <- PIRLS(d, ex, w, nx, bdeg, 2, 1, lambda1, lambda2)
      X.MM <- bspline(x, nx, bdeg)$BB
      for.t2s1 <- c(X.MM %*% theta.MM$theta)[(tfit+1):ttot]
      bic.t2s1[i] <- theta.MM$Bic
      rmse.t2s1[i] <- (sqrt((1/n.for)*(sum((exp(for.t2s1)-r.obs)^2))))
      mape.t2s1[i] <- ((100/n.for)*sum(abs((r.obs-exp(for.t2s1))/r.obs)))
      
      print(i)
    }
    
    country.db[[j]] <- 
      cbind(
        na.omit(bic.t1s1),na.omit(bic.t2s2),na.omit(bic.t1s2),na.omit(bic.t2s1),
        na.omit(rmse.t1s1),na.omit(rmse.t2s2),na.omit(rmse.t1s2),na.omit(rmse.t2s1),
        na.omit(mape.t1s1),na.omit(mape.t2s2),na.omit(mape.t1s2),na.omit(mape.t2s1))
    print(country.names[j])
    
  }
  
  return(country.db)
  
}


#####-------------- Accuracy - 5 years rolling window --------------#####

# Starting month of the rolling window
window <- seq(1,229,12)
# Start and end month of the forecasts in each window
tfit <- 60
ttot <- 72
# Weight vector (1 estimation set, 0 forecast set)
w <- rep(1,ttot)
w[(tfit+1):ttot] <- 0
# Degree of basis and number of basis (take 6*2 basis)
bdeg = 3
nx = 12

# Accuracy measures
db <- as.data.frame(matrix(NA,length(window),12))
country.db <- replicate(25, db, simplify=F)
fit.db <- accuracy(rates, country.names, country.db)
db.5y.av <- rbind(
  colMeans(fit.db[[1]]), colMeans(fit.db[[2]]), colMeans(fit.db[[3]]), 
  colMeans(fit.db[[4]]), colMeans(fit.db[[5]]), colMeans(fit.db[[6]]),
  colMeans(fit.db[[7]]), colMeans(fit.db[[8]]), colMeans(fit.db[[9]]), 
  colMeans(fit.db[[10]]), colMeans(fit.db[[11]]), colMeans(fit.db[[12]]),
  colMeans(fit.db[[13]]), colMeans(fit.db[[14]]), colMeans(fit.db[[15]]), 
  colMeans(fit.db[[16]]), colMeans(fit.db[[17]]), colMeans(fit.db[[18]]),
  colMeans(fit.db[[19]]), colMeans(fit.db[[20]]), colMeans(fit.db[[21]]), 
  colMeans(fit.db[[22]]), colMeans(fit.db[[23]]), colMeans(fit.db[[24]]),
  colMeans(fit.db[[25]]))
colnames(db.5y.av) <- accuracy.names
rownames(db.5y.av) <- country.names


#####-------------- Accuracy - 10 years rolling window --------------#####

# Starting month of the rolling window
window <- seq(1,169,12)
# Start and end month of the forecasts in each window
tfit <- 120
ttot <- 132
# Weight vector (1 estimation set, 0 forecast set)
w <- rep(1,ttot)
w[(tfit+1):ttot] <- 0
# Degree of basis and number of basis (take 11*2 basis)
bdeg = 3
nx = 22

# Accuracy measures
db <- as.data.frame(matrix(NA,length(window),12))
country.db <- replicate(25, db, simplify=F)
fit.db <- accuracy(rates, country.names, country.db)
db.10y.av <- rbind(
  colMeans(fit.db[[1]]), colMeans(fit.db[[2]]), colMeans(fit.db[[3]]), 
  colMeans(fit.db[[4]]), colMeans(fit.db[[5]]), colMeans(fit.db[[6]]),
  colMeans(fit.db[[7]]), colMeans(fit.db[[8]]), colMeans(fit.db[[9]]), 
  colMeans(fit.db[[10]]), colMeans(fit.db[[11]]), colMeans(fit.db[[12]]),
  colMeans(fit.db[[13]]), colMeans(fit.db[[14]]), colMeans(fit.db[[15]]), 
  colMeans(fit.db[[16]]), colMeans(fit.db[[17]]), colMeans(fit.db[[18]]),
  colMeans(fit.db[[19]]), colMeans(fit.db[[20]]), colMeans(fit.db[[21]]), 
  colMeans(fit.db[[22]]), colMeans(fit.db[[23]]), colMeans(fit.db[[24]]),
  colMeans(fit.db[[25]]))
colnames(db.10y.av) <- accuracy.names
rownames(db.10y.av) <- country.names


#####-------------- Save environment death rates --------------#####

# save(fit.db, db.5y.av, file="fig/penalty.EU.R5y.RData")
# save(fit.db, db.10y.av, file="fig/penalty.EU.R10y.RData")


######################################################
##########        SAVE LATEX TABLES         ##########
######################################################
#####-------------- Death rates --------------#####

load("fig/accuracy/penalty.EU.R5y.RData")
load("fig/accuracy/penalty.EU.R10y.RData")

# Remember to change [t] in [c]
# Change textsuperscript{a} in textsuperscript{*}

rmse.av <- cbind(db.5y.av[,5:12],db.10y.av[,5:12])
rmse.av <- rmse.av[ order(row.names(rmse.av)), ]
rmse.av[,c(1:4,9:12)] <- round(rmse.av[,c(1:4,9:12)]*1000,2)
rmse.av[,c(5:8,13:16)] <- round(rmse.av[,c(5:8,13:16)],2)
rmse.av <- as.data.frame(rmse.av) 
Country = rownames(rmse.av)
rownames(rmse.av) <- c(1:25)
rmse.av <- cbind(Country, rmse.av) 
colnames(rmse.av) <- c("Country",rep(c("(1)","(2)","(3)","(4)"),4))
rmse.av <- rbind(rmse.av,c("No. *",15,0,0,10,12,0,0,13,6,0,0,19,3,1,0,21))
caption= "Effect of choice of the order of the penalty on the forecasts. Mean RMSE and MAPE on CDRs in 25 European countries for multiple fitting periods based on a rolling-window scheme. Four combinations of penalties  are considered: (1) order 1 on trend and seasonality, (2) order 2 on the trend and seasonality, (3) order 1 on the trend and 2 on the seasonality, (4) order 2 on the trend and order 1 on the seasonality."
# for(i in c(2:5,10:13)){ rmse.av[,i] <- comma(rmse.av[,i]) }
rmse.av %>%
  kbl(align="lcccccccccc", longtable=T, booktabs=T,
      caption=caption, escape=F, format="latex", label="best2") %>%
  kable_styling(font_size = 7.5) %>%
  add_header_above(c(" "=1,"RMSE"=4,"MAPE"=4,"RMSE"=4,"MAPE"=4)) %>%
  add_header_above(c(" "=1,"5 years series"=8,"10 years series"=8)) %>%
  column_spec(1, border_right=T) %>% column_spec(5, border_right=T) %>%
  column_spec(9, border_right=T) %>% column_spec(13, border_right=T) %>%
  footnote(alphabet = c("Number of countries in which the model performs the best"))


######################################################
##########    PLOT EFFECT OF THE PENALTY    ##########
######################################################
#####-------------- Directory and packages --------------#####

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
load("fig/data/series.EU.RData")

library(tidyverse)
require(ggplot2)
library(splines)
library(zoo)
library(data.table)
library(readxl)
library(foreign)
library(magic)

# Load function for forecasting
source("fun/ForFun total.R")

#####-------------- Parameters choice --------------#####

swe1 <- rates %>% filter(country=="Sweden", date>="2010-01-01", 
                         date<="2022-12-01")
d <- (swe1 %>% arrange(date))$deaths
ex <- (swe1 %>% arrange(date))$exposure
x <- 1:length(d)
period <- sort(unique(swe1$date))
tfit <- 121 # ending month of the fit
w <- rep(1,length(period)) # weight vector (1 estimation set)
w[(tfit+1):length(period)] <- 0 # weight vector (0 forecast set)
d.fit <- d[w==1] # subset estimation dataset
ex.fit <- ex[w==1]
nx = 24 # number of regions, 12years (10fit,2forecast) * 2basis
bdeg = 3 # degree of basis


#####-------------- Trend --------------#####

##### Estimate the modulation model for Spain (pord1=1, pord2=1)

# Grid search for lambda1 and lambda2 (MM) and lambda (SM)
lambda.seq <- 10^seq(5,7,length.out=30)
mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
mape.s.MM <- rep(NA,300)
for (h in 1:length(lambda.seq)){
  for(g in 1:length(lambda.seq)){
    for(k in seq(13,109,12)){ 
      d.sel <- d.fit[1:(k+11)]
      ex.sel <- ex.fit[1:(k+11)]
      w.sel <- rep(1,k+11)
      w.sel[k:(k+11)] <- 0
      nx.sel <- (floor(k/12)+1)*2
      Output.MM <- PIRLS(d.sel, ex.sel, w.sel, nx.sel, bdeg, 1, 1, 
                         lambda.seq[h], lambda.seq[g])$theta
      x.s <- 1:length(d.sel)
      X.s <- bspline(x.s, nx.sel, bdeg)
      eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
      r.obs.s <- (d.fit/ex.fit)[k:(k+11)]
      mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
    }
    mape.mat[h,g] <- mean(na.omit(mape.s.MM))
  }}
mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
lambda1 = lambda.seq[mape.min.MM[1]]
lambda2 = lambda.seq[mape.min.MM[2]]

theta.MM <- PIRLS(d, ex, w, nx, bdeg, 1, 1, lambda1, lambda2)
B.MM <- bspline(x, nx, bdeg)$BB
vt.MM <- c(B.MM[,1:(nx+bdeg)] %*% theta.MM$theta[1:(nx+bdeg)])
db.MM1 <- data.frame(
  pred=vt.MM,eta.up=NA,eta.lo=NA,
  year=period,year.fit=c(period[w==1],rep(NA,length(w[w==0]))),
  rates=d/ex,rates.fit=c(d.fit[w==1]/ex.fit[w==1],rep(NA,length(w[w==0]))),
  modify="trend",model="pord1=1, pord2=1")

##### Estimate the modulation model for Spain (pord1=2, pord2=1)

# Grid search for lambda1 and lambda2 (MM) and lambda (SM)
lambda.seq <- 10^seq(5,7,length.out=30)
mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
mape.s.MM <- rep(NA,300)
for (h in 1:length(lambda.seq)){
  for(g in 1:length(lambda.seq)){
    for(k in seq(13,109,12)){ 
      d.sel <- d.fit[1:(k+11)]
      ex.sel <- ex.fit[1:(k+11)]
      w.sel <- rep(1,k+11)
      w.sel[k:(k+11)] <- 0
      nx.sel <- (floor(k/12)+1)*2
      Output.MM <- PIRLS(d.sel, ex.sel, w.sel, nx.sel, bdeg, 2, 1, 
                         lambda.seq[h], lambda.seq[g])$theta
      x.s <- 1:length(d.sel)
      X.s <- bspline(x.s, nx.sel, bdeg)
      eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
      r.obs.s <- (d.fit/ex.fit)[k:(k+11)]
      mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
    }
    mape.mat[h,g] <- mean(na.omit(mape.s.MM))
  }}
mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
lambda1 = lambda.seq[mape.min.MM[1]]
lambda2 = lambda.seq[mape.min.MM[2]]

theta.MM <- PIRLS(d, ex, w, nx, bdeg, 2, 1, lambda1, lambda2)
B.MM <- bspline(x, nx, bdeg)$BB
vt.MM <- c(B.MM[,1:(nx+bdeg)] %*% theta.MM$theta[1:(nx+bdeg)])
db.MM2 <- data.frame(
  pred=vt.MM,eta.up=NA,eta.lo=NA,
  year=period,year.fit=c(period[w==1],rep(NA,length(w[w==0]))),
  rates=d/ex,rates.fit=c(d.fit[w==1]/ex.fit[w==1],rep(NA,length(w[w==0]))),
  modify="trend",model="pord1=2, pord2=1")

##### Estimate the modulation model for Spain (pord1=3, pord2=1)

# Grid search for lambda1 and lambda2 (MM) and lambda (SM)
lambda.seq <- 10^seq(5,7,length.out=30)
mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
mape.s.MM <- rep(NA,300)
for (h in 1:length(lambda.seq)){
  for(g in 1:length(lambda.seq)){
    for(k in seq(13,109,12)){ 
      d.sel <- d.fit[1:(k+11)]
      ex.sel <- ex.fit[1:(k+11)]
      w.sel <- rep(1,k+11)
      w.sel[k:(k+11)] <- 0
      nx.sel <- (floor(k/12)+1)*2
      Output.MM <- PIRLS(d.sel, ex.sel, w.sel, nx.sel, bdeg, 3, 1, 
                         lambda.seq[h], lambda.seq[g])$theta
      x.s <- 1:length(d.sel)
      X.s <- bspline(x.s, nx.sel, bdeg)
      eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
      r.obs.s <- (d.fit/ex.fit)[k:(k+11)]
      mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
    }
    mape.mat[h,g] <- mean(na.omit(mape.s.MM))
  }}
mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
lambda1 = lambda.seq[mape.min.MM[1]]
lambda2 = lambda.seq[mape.min.MM[2]]

theta.MM <- PIRLS(d, ex, w, nx, bdeg, 3, 1, lambda1, lambda2)
B.MM <- bspline(x, nx, bdeg)$BB
vt.MM <- c(B.MM[,1:(nx+bdeg)] %*% theta.MM$theta[1:(nx+bdeg)])
db.MM3 <- data.frame(
  pred=vt.MM,eta.up=NA,eta.lo=NA,
  year=period,year.fit=c(period[w==1],rep(NA,length(w[w==0]))),
  rates=d/ex,rates.fit=c(d.fit[w==1]/ex.fit[w==1],rep(NA,length(w[w==0]))),
  modify="trend",model="pord1=3, pord2=1")

#####-------------- Seasonality --------------#####

##### Estimate the modulation model for Spain (pord1=2, pord2=1)

# Grid search for lambda1 and lambda2 (MM) and lambda (SM)
lambda.seq <- 10^seq(5,7,length.out=30)
mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
mape.s.MM <- rep(NA,300)
for (h in 1:length(lambda.seq)){
  for(g in 1:length(lambda.seq)){
    for(k in seq(13,109,12)){ 
      d.sel <- d.fit[1:(k+11)]
      ex.sel <- ex.fit[1:(k+11)]
      w.sel <- rep(1,k+11)
      w.sel[k:(k+11)] <- 0
      nx.sel <- (floor(k/12)+1)*2
      Output.MM <- PIRLS(d.sel, ex.sel, w.sel, nx.sel, bdeg, 2, 1, 
                         lambda.seq[h], lambda.seq[g])$theta
      x.s <- 1:length(d.sel)
      X.s <- bspline(x.s, nx.sel, bdeg)
      eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
      r.obs.s <- (d.fit/ex.fit)[k:(k+11)]
      mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
    }
    mape.mat[h,g] <- mean(na.omit(mape.s.MM))
  }}
mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
lambda1 = lambda.seq[mape.min.MM[1]]
lambda2 = lambda.seq[mape.min.MM[2]]

theta.MM <- PIRLS(d, ex, w, nx, bdeg, 2, 1, lambda1, lambda2)
B.MM <- bspline(x, nx, bdeg)$BB
eta.MM <- c(B.MM %*% theta.MM$theta)
c.val <- qnorm(0.975)
eta.up.MM <- eta.MM + c.val * theta.MM$se.eta
eta.lo.MM <- eta.MM - c.val * theta.MM$se.eta
db.MM4 <- data.frame(
  pred=eta.MM,eta.up=eta.up.MM,eta.lo=eta.lo.MM,
  year=period,year.fit=c(period[w==1],rep(NA,length(w[w==0]))),
  rates=d/ex,rates.fit=c(d.fit[w==1]/ex.fit[w==1],rep(NA,length(w[w==0]))),
  modify="season",model="pord1=2, pord2=1")

##### Estimate the modulation model for Spain (pord1=2, pord2=1)

# Grid search for lambda1 and lambda2 (MM) and lambda (SM)
lambda.seq <- 10^seq(5,7,length.out=30)
mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
mape.s.MM <- rep(NA,300)
for (h in 1:length(lambda.seq)){
  for(g in 1:length(lambda.seq)){
    for(k in seq(13,109,12)){ 
      d.sel <- d.fit[1:(k+11)]
      ex.sel <- ex.fit[1:(k+11)]
      w.sel <- rep(1,k+11)
      w.sel[k:(k+11)] <- 0
      nx.sel <- (floor(k/12)+1)*2
      Output.MM <- PIRLS(d.sel, ex.sel, w.sel, nx.sel, bdeg, 2, 2, 
                         lambda.seq[h], lambda.seq[g])$theta
      x.s <- 1:length(d.sel)
      X.s <- bspline(x.s, nx.sel, bdeg)
      eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
      r.obs.s <- (d.fit/ex.fit)[k:(k+11)]
      mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
    }
    mape.mat[h,g] <- mean(na.omit(mape.s.MM))
  }}
mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
lambda1 = lambda.seq[mape.min.MM[1]]
lambda2 = lambda.seq[mape.min.MM[2]]

theta.MM <- PIRLS(d, ex, w, nx, bdeg, 2, 2, lambda1, lambda2)
B.MM <- bspline(x, nx, bdeg)$BB
eta.MM <- c(B.MM %*% theta.MM$theta)
c.val <- qnorm(0.975)
eta.up.MM <- eta.MM + c.val * theta.MM$se.eta
eta.lo.MM <- eta.MM - c.val * theta.MM$se.eta
db.MM5 <- data.frame(
  pred=eta.MM,eta.up=eta.up.MM,eta.lo=eta.lo.MM,
  year=period,year.fit=c(period[w==1],rep(NA,length(w[w==0]))),
  rates=d/ex,rates.fit=c(d.fit[w==1]/ex.fit[w==1],rep(NA,length(w[w==0]))),
  modify="season",model="pord1=2, pord2=2")

##### Estimate the modulation model for Spain (pord1=3, pord2=1)

# Grid search for lambda1 and lambda2 (MM) and lambda (SM)
lambda.seq <- 10^seq(5,7,length.out=30)
mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
mape.s.MM <- rep(NA,300)
for (h in 1:length(lambda.seq)){
  for(g in 1:length(lambda.seq)){
    for(k in seq(13,109,12)){ 
      d.sel <- d.fit[1:(k+11)]
      ex.sel <- ex.fit[1:(k+11)]
      w.sel <- rep(1,k+11)
      w.sel[k:(k+11)] <- 0
      nx.sel <- (floor(k/12)+1)*2
      Output.MM <- PIRLS(d.sel, ex.sel, w.sel, nx.sel, bdeg, 2, 3, 
                         lambda.seq[h], lambda.seq[g])$theta
      x.s <- 1:length(d.sel)
      X.s <- bspline(x.s, nx.sel, bdeg)
      eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
      r.obs.s <- (d.fit/ex.fit)[k:(k+11)]
      mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
    }
    mape.mat[h,g] <- mean(na.omit(mape.s.MM))
  }}
mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
lambda1 = lambda.seq[mape.min.MM[1]]
lambda2 = lambda.seq[mape.min.MM[2]]

theta.MM <- PIRLS(d, ex, w, nx, bdeg, 2, 3, lambda1, lambda2)
B.MM <- bspline(x, nx, bdeg)$BB
eta.MM <- c(B.MM %*% theta.MM$theta)
c.val <- qnorm(0.975)
eta.up.MM <- eta.MM + c.val * theta.MM$se.eta
eta.lo.MM <- eta.MM - c.val * theta.MM$se.eta
db.MM6 <- data.frame(
  pred=eta.MM,eta.up=eta.up.MM,eta.lo=eta.lo.MM,
  year=period,year.fit=c(period[w==1],rep(NA,length(w[w==0]))),
  rates=d/ex,rates.fit=c(d.fit[w==1]/ex.fit[w==1],rep(NA,length(w[w==0]))),
  modify="season",model="pord1=2, pord2=3")

#####-------------- Plot --------------#####

data.fig <- rbind(db.MM1,db.MM2,db.MM3, db.MM4,db.MM5,db.MM6) %>%
  mutate(
    type = case_when(
      modify == "trend" & model == "pord1=1, pord2=1" ~ "pord=1",
      modify == "trend" & model == "pord1=2, pord2=1" ~ "pord=2",
      modify == "trend" & model == "pord1=3, pord2=1" ~ "pord=3",
      modify == "season" & model == "pord1=2, pord2=1" ~ "pord=1",
      modify == "season" & model == "pord1=2, pord2=2" ~ "pord=2",
      modify == "season" & model == "pord1=2, pord2=3" ~ "pord=3")) %>%
  mutate(modify = factor(modify, 
                         levels=c("trend","season"), 
                         labels=c("trend"="trend","season"="seasonality")))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 3
cols = gg_color_hue(n)
codes.col = c("pord=1"=cols[1], "pord=2"=cols[2], "pord=3"=cols[3])

fig <- data.fig %>%
  filter(modify == "trend") %>%
  ggplot() +
  geom_line(aes(x=year, y=log(rates)), color="red", size=0.5) +  
  geom_line(aes(x=year.fit, y=log(rates.fit)), color="grey", size=0.5) +
  geom_line(aes(x=year, y=pred, group=type, colour=type), linetype="dashed", size=0.5) +
  # geom_line(aes(x=year, y=pred, group=type, colour=type), size=1) +
  # geom_smooth(data=(data.fig %>% filter(modify=="seasonality")),
  #             aes(x=year, y=pred, ymin=unname(eta.lo), ymax=unname(eta.up),
  #                 fill=type, colour=type, group=type),
  #             linetype="solid", stat='identity') +
  # facet_grid(vars(modify), scales = "free") +
  scale_x_date(date_labels="%Y", date_breaks="2 year",
               expand=c(0.05,0.05)) +  
  scale_colour_manual(values=codes.col) +
  labs(title="Sweden", x="Month", y="Log death rates") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 9),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.position="bottom",
    legend.key.width=unit(2.5, "cm"),
    plot.title = element_text(hjust = 0.5, size = 9))

ggsave("fig/accuracy/FigEU 77.png", plot=fig, width=7, height=4, bg="white")

