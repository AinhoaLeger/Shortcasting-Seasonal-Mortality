
######################################################
##########        ACCURACY FOR RATES        ##########
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
library(kableExtra)
library(knitr)

# Load function for forecasting
source("fun/ForFun total.R")


#####-------------- Function for rates --------------#####

period <- sort(unique(rates$date))

# Names of the sub-population and accuracy measures
accuracy.names <- c("bic.SP","bic.MM","bic.SM",
                    "rmse.SP","rmse.MM","rmse.SM",
                    "mape.SP","mape.MM","mape.SM")
country.names= c("Bulgaria","Czechia","Denmark","Germany",
                 "Estonia","Ireland","Greece","Spain","France",
                 "Croatia","Italy","Lithuania","Luxembourg",
                 "Hungary","Netherlands","Austria","Poland",
                 "Portugal","Romania","Slovenia","Finland",
                 "Sweden","Iceland","Norway","Switzerland")

accuracy <- function(data, country.names, country.db) {
  
  bic.SP <- bic.MM <- bic.SM <- rmse.SP <- rmse.MM <- rmse.SM <-
    mape.SP <- mape.MM <- mape.SM <-  c()
  
  for(j in 1:length(country.names)){
    for(i in window){ 
      
      d <- (data %>% filter(country==country.names[j]))$deaths[i:(i+ttot-1)]
      ex <- (data %>% filter(country==country.names[j]))$exposure[i:(i+ttot-1)]
      x <- 1:length(d)
      d.fit <- d[w==1] # subset estimation dataset
      ex.fit <- ex[w==1]
      n.for <- 12 # we produce 12 forecast in a year
      r.obs <- (d/ex)[(tfit+1):ttot] # subset prediction dataset
      
      #####----------------- Poisson Serfling -----------------######
      
      df = data.frame(x=x, d=d, cos=cos(((2*pi)/12)*x), sin=sin(((2*pi)/12)*x))
      XS <- model.matrix(d ~ x + cos + sin, data = df)
      Dinit <- matrix(log( d+1 ), nrow = nrow(XS), ncol = 1)
      betainit <- solve( t(XS) %*% XS, t(XS) %*% Dinit)
      betaS <- IRLS(XS, w, d, ex, betainit)
      eta.SP <- c(betaS$theta[1] + betaS$theta[2]*x + 
        betaS$theta[3]*cos(((2*pi)/12)*x) + 
          betaS$theta[4]*sin(((2*pi)/12)*x))[(tfit+1):ttot]
      bic.SP[i] <- betaS$Bic
      rmse.SP[i] <- (sqrt((1/n.for)*(sum((exp(eta.SP)-r.obs)^2)))) 
      mape.SP[i] <- ((100/n.for)*sum(abs((r.obs-exp(eta.SP))/r.obs))) 
      
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
          for(k in c(25,37,49)){ # 10 years: 73,85,97,109
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

      theta.MM <- PIRLS(d, ex, w, nx, bdeg, 2, 1, lambda1, lambda2)
      X.MM <- bspline(x, nx, bdeg)$BB
      eta.MM <- c(X.MM %*% theta.MM$theta)[(tfit+1):ttot]
      bic.MM[i] <- theta.MM$Bic
      rmse.MM[i] <- (sqrt((1/n.for)*(sum((exp(eta.MM)-r.obs)^2))))
      mape.MM[i] <- ((100/n.for)*sum(abs((r.obs-exp(eta.MM))/r.obs)))
      
      #####----------------- Smooth trend model -----------------######
      
      theta.SM <- PIRLS2(d, ex, w, nx, bdeg, 2, lambda)
      X.SM <- bspline(x, nx, bdeg)$BB2
      eta.SM <- c(X.SM %*% theta.SM$theta)[(tfit+1):ttot]      
      bic.SM[i] <- theta.SM$Bic
      rmse.SM[i] <- (sqrt((1/n.for)*(sum((exp(eta.SM)-r.obs)^2))))
      mape.SM[i] <- ((100/n.for)*sum(abs((r.obs-exp(eta.SM))/r.obs)))
      
      rm(betaS,eta.SP, Output.MM,theta.MM,eta.MM, Output.SM,theta.SM,eta.SM)
      print(i)
    }
    
    country.db[[j]] <- 
      cbind(na.omit(bic.SP),na.omit(bic.MM),na.omit(bic.SM),
            na.omit(rmse.SP),na.omit(rmse.MM),na.omit(rmse.SM),
            na.omit(mape.SP),na.omit(mape.MM),na.omit(mape.SM))
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
db <- as.data.frame(matrix(NA,length(window),9))
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
db <- as.data.frame(matrix(NA,length(window),9))
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

# save(fit.db, db.5y.av, file="fig/accuracy/accuracy.EU.R5y.RData")
# save(fit.db, db.10y.av, file="fig/accuracy/accuracy.EU.R10y.RData")

######################################################
##########        SAVE LATEX TABLES         ##########
######################################################
#####-------------- Death rates --------------#####

load("fig/accuracy/accuracy.EU.R5y.RData")
load("fig/accuracy/accuracy.EU.R10y.RData")

# Remember to change [t] in [c]
# Change textsuperscript{a} in textsuperscript{*}

Bic.av <- cbind(db.5y.av[,1:3],db.10y.av[,1:3])
Bic.av <- Bic.av[ order(row.names(Bic.av)), ]
Bic.av <- round(Bic.av)
Bic.av <- as.data.frame(Bic.av) 
Country = rownames(Bic.av)
rownames(Bic.av) <- c(1:25)
Bic.av <- cbind(Country, Bic.av) 
colnames(Bic.av) <- c("Country",rep(c("SP","SP-STSS","SP-STFS"),2))
Bic.av <- rbind(Bic.av,c("No. countries *",4,20,1,4,21,0))
caption= "Mean BIC on CDRs in 25 European countries for multiple fitting periods based on a rolling-window scheme. Three models (SP, SP-STSS and SP-STFS) and two lengths of the time series are compared (5 and 10 years)."
for(i in 2:7){ Bic.av[,i] <- comma(Bic.av[,i]) }
Bic.av %>%
  kbl(align="lcccccccccc", longtable=T, booktabs=T,
      caption=caption, escape=F, format="latex", label="bic") %>%
  kable_styling(font_size = 8) %>%
  add_header_above(c(" "=1,"5 years series"=3,"10 years series"=3)) %>%
  column_spec(1, border_right=T) %>% column_spec(4, border_right=T) %>%
  footnote(alphabet = c("Number of countries in which the model performs the best"))

# Remember to change [t] in [c]
# Change textsuperscript{a} in textsuperscript{*}

rmse.av <- cbind(db.5y.av[,4:9],db.10y.av[,4:9])
rmse.av <- rmse.av[ order(row.names(rmse.av)), ]
rmse.av[,c(1:3,7:9)] <- round(rmse.av[,c(1:3,7:9)]*1000,2)
rmse.av[,c(4:6,10:12)] <- round(rmse.av[,c(4:6,10:12)],2)
rmse.av <- as.data.frame(rmse.av)
Country = rownames(rmse.av)
rownames(rmse.av) <- c(1:25)
rmse.av <- cbind(Country, rmse.av)
colnames(rmse.av) <- c("Country",rep(c("SP","STSS","STFS"),4))
rmse.av <- rbind(rmse.av,c("No. countries *",17,2,6,18,1,6,7,4,14,9,3,13))
caption= "Mean RMSE and MAPE on one-year ahead CDRs (x 1000) in 25 European countries for multiple fitting periods based on a rolling-window scheme. Three models (SP, SP-STSS and SP-STFS) and two lengths of the time series are compared (5 and 10 years)."
rmse.av %>%
  kbl(align="lcccccccccc", longtable=T, booktabs=T,
      caption=caption, escape=F, format="latex", label="mape") %>%
  kable_styling(font_size = 8) %>%
  add_header_above(c(" "=1,"RMSE"=3,"MAPE"=3,"RMSE"=3,"MAPE"=3)) %>%
  add_header_above(c(" "=1,"5 years series"=6,"10 years series"=6)) %>%
  column_spec(1, border_right=T) %>% column_spec(4, border_right=T) %>%
  column_spec(7, border_right=T) %>% column_spec(10, border_right=T) %>%
  footnote(alphabet = c("Number of countries in which the model performs the best"))

