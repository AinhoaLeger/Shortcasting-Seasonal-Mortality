
######################################################
##########            PARAMETERS            ##########
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
                         date<="2020-02-01")
d <- (swe1 %>% arrange(date))$deaths
e <- (swe1 %>% arrange(date))$exposure
x <- 1:length(d)
period <- sort(unique(swe1$date))
w <- rep(1,length(period))
off <- log(e)
nx = 20 # number of regions (10 years * 2 regions)
bdeg = 3 # degree of basis

##### Grid search for lambda

lambda.seq <- 10^seq(4,7,length.out=30)
mape.vec <- rep(NA, length(lambda.seq)) 
mape.s.SM <- rep(NA,300)
for(g in 1:length(lambda.seq)){
  # Rolling window on the estimation set
  for(k in seq(13,109,12)){
    d.sel <- d[1:(k+11)]
    e.sel <- e[1:(k+11)]
    w.sel <- rep(1,k+11)
    w.sel[k:(k+11)] <- 0
    nx.sel <- (floor(k/12)+1)*2
    Output.SM <- PIRLS2(d.sel, e.sel, w.sel, nx.sel, bdeg, 2,
                        lambda.seq[g])$theta
    x.s <- 1:length(d.sel)
    X.s <- bspline(x.s, nx.sel, bdeg)
    eta.s.SM <- c(X.s$BB2 %*% Output.SM)[k:(k+11)]
    r.obs.s <- (d/e)[k:(k+11)]
    mape.s.SM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.SM))/r.obs.s)))
  }
  mape.vec[g] <- mean(na.omit(mape.s.SM))
}
mape.min.SM <- which(mape.vec == min(mape.vec), arr.ind = TRUE)
lambda = lambda.seq[mape.min.SM[1]]

##### Estimate the smooth trend model

theta.SM <- PIRLS2(d, e, w, nx, bdeg, 2, lambda)
B <- bspline(1:length(period), nx, bdeg)$B
B.SM <- bspline(x, nx, bdeg)$BB2
# BB2: 51 trend coefficients (48+3), 2 seasonality coefficients
vt.sca <- t(t(B.SM[,1:(nx+bdeg)]) * as.vector(theta.SM$theta[1:(nx+bdeg)]))
vt.add <- c(B.SM[,1:(nx+bdeg)] %*% theta.SM$theta[1:(nx+bdeg)])


##### Plots basis

png(file = "fig/forecast/FigEU 14.png", width=800, height=700) 

par(mfrow=c(2,1), oma = c(2,2,0,0) + 1, mar = c(0,3.7,1,1) + 1)
# Plot data and B-splines basis
matplot(period, B, type="l", lwd=2, cex.axis=1.5, ylim=c(0,0.7), ylab="")
# Plot data, rescaled and added basis B-splines basis, and fit model
# matplot(period, log(d/e), type="l", cex.axis=2, ylim=c(-3.3,0),
#         ylab="Death rates (log)", cex.lab=2.5)
# matplot(period, vt.sca, type="l", lwd=2, cex.axis=2, ylab="", add=TRUE)
# Plot data, rescaled and added basis B-splines basis, and fit model
matplot(period, log(d/e), type="l", cex.axis=1.5,
        ylab="Death rates (log)", cex.lab=2)
matplot(period, vt.add, type="l", lwd=2, col="cyan3", cex.axis=2, add=TRUE)

dev.off()


######################################################
##########         MODELLING RATES          ##########
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

# Load function for forecasting
source("fun/ForFun total.R")


#####-------------- Fit and forecast --------------#####

# Select 10 year for fitting (common period)
eu <- rates %>% filter(date >= "2010-01-01", date <= "2022-06-01")
period <- sort(unique(eu$date))
country.list=unique(eu$country)
# Start and end month of the forecasts in each window
tfit <- 121
# Weight vector (1 estimation set, 0 forecast set)
w <- rep(1,length(period))
w[(tfit+1):length(period)] <- 0
# Number and degree of basis = 12 years (10 fit, 2 forecast) * 2 basis
nx = 24
bdeg = 3


db.eu.SP <- db.eu.MM <- db.eu.SM <- list()

for(i in 1:length(country.list)){
  
  d <- (eu %>% filter(country == country.list[i]) %>% arrange(date))$deaths
  e <- (eu %>% filter(country == country.list[i]) %>% arrange(date))$exposure
  x <- 1:length(d)
  d.fit <- d[w==1] # subset estimation dataset
  e.fit <- e[w==1]
  
  #####----------------- Poisson Serfling -----------------######
  
  df = data.frame(x=x, d=d, cos=cos(((2*pi)/12)*x), sin=sin(((2*pi)/12)*x))
  XT <- model.matrix(d ~ x, data = df)
  XS <- model.matrix(d ~ x + cos + sin, data = df)
  Dinit <- matrix(log( d+1 ), nrow = nrow(XS), ncol = 1)
  betainitT <- solve( t(XT) %*% XT, t(XT) %*% Dinit)
  betainitS <- solve( t(XS) %*% XS, t(XS) %*% Dinit)
  betaT <- IRLS(XT, w, d, e, betainitT)
  betaS <- IRLS(XS, w, d, e, betainitS)
  vt.SP <- betaT$theta[1] + betaT$theta[2]*x
  eta.SP <- c(betaS$theta[1] + betaS$theta[2]*x + 
                betaS$theta[3]*cos(((2*pi)/12)*x) + 
                betaS$theta[4]*sin(((2*pi)/12)*x))
  c.val <- qnorm(0.975)
  eta.up.SP <- eta.SP + c.val * betaS$se.eta
  eta.lo.SP <- eta.SP - c.val * betaS$se.eta
  
  db.eu.SP[[i]] <- data.frame(vt=vt.SP, eta=eta.SP, 
                              eta.up=eta.up.SP, eta.lo=eta.lo.SP)
  db.eu.SP[[i]]$country <- country.list[i]
  db.eu.SP[[i]]$rates = d/e
  db.eu.SP[[i]]$rates.fit = c(d[w!=0]/e[w!=0], rep(NA,length(w[w==0])))
  db.eu.SP[[i]]$year <- period
  db.eu.SP[[i]]$year.fit = c(period[w!=0], rep(NA,length(w[w==0])))
  db.eu.SP[[i]]$model = "Poisson Serfling"
  
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
      for(k in seq(49,109,12)){ 
        d.sel <- d.fit[1:(k+11)]
        e.sel <- e.fit[1:(k+11)]
        w.sel <- rep(1,k+11)
        w.sel[k:(k+11)] <- 0
        nx.sel <- (floor(k/12)+1)*2
        Output.MM <- PIRLS(d.sel, e.sel, w.sel, nx.sel, bdeg, 2, 1, 
                           lambda.seq[h], lambda.seq[g])$theta
        Output.SM <- PIRLS2(d.sel, e.sel, w.sel, nx.sel, bdeg, 2,
                            lambda.seq[g])$theta
        x.s <- 1:length(d.sel)
        X.s <- bspline(x.s, nx.sel, bdeg)
        eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
        eta.s.SM <- c(X.s$BB2 %*% Output.SM)[k:(k+11)]
        r.obs.s <- (d.fit/e.fit)[k:(k+11)]
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
  
  theta.MM <- PIRLS(d, e, w, nx, bdeg, 2, 1, lambda1, lambda2)
  X.MM <- bspline(x, nx, bdeg)$BB
  vt.MM <- c(X.MM[,1:(nx+bdeg)] %*% theta.MM$theta[1:(nx+bdeg)])
  eta.MM <- c(X.MM %*% theta.MM$theta)
  c.val <- qnorm(0.975)
  eta.up.MM <- eta.MM + c.val * theta.MM$se.eta
  eta.lo.MM <- eta.MM - c.val * theta.MM$se.eta
  
  db.eu.MM[[i]] <- data.frame(vt=vt.MM, eta=eta.MM, 
                              eta.up=eta.up.MM, eta.lo=eta.lo.MM)
  db.eu.MM[[i]]$country <- country.list[i]
  db.eu.MM[[i]]$rates = d/e
  db.eu.MM[[i]]$rates.fit = c(d[w!=0]/e[w!=0], rep(NA,length(w[w==0])))
  db.eu.MM[[i]]$year <- period
  db.eu.MM[[i]]$year.fit = c(period[w!=0], rep(NA,length(w[w==0])))
  db.eu.MM[[i]]$model = "Modulation model"
  
  #####----------------- Smooth trend model -----------------######
  
  theta.SM <- PIRLS2(d, e, w, nx, bdeg, 2, lambda)
  X.SM <- bspline(x, nx, bdeg)$BB2
  vt.SM <- c(X.SM[,1:(nx+bdeg)] %*% theta.SM$theta[1:(nx+bdeg)])
  eta.SM <- c(X.SM %*% theta.SM$theta)
  c.val <- qnorm(0.975)
  eta.up.SM <- eta.SM + c.val * theta.SM$se.eta
  eta.lo.SM <- eta.SM - c.val * theta.SM$se.eta
  
  db.eu.SM[[i]] <- data.frame(vt=vt.SM, eta=eta.SM, 
                              eta.up=eta.up.SM, eta.lo=eta.lo.SM)
  db.eu.SM[[i]]$country <- country.list[i]
  db.eu.SM[[i]]$rates = d/e
  db.eu.SM[[i]]$rates.fit = c(d[w!=0]/e[w!=0], rep(NA,length(w[w==0])))
  db.eu.SM[[i]]$year <- period
  db.eu.SM[[i]]$year.fit = c(period[w!=0], rep(NA,length(w[w==0])))
  db.eu.SM[[i]]$model = "Smooth trend model"
  
  rm(betaS,eta.SP, Output.MM,theta.MM,eta.MM, Output.SM,theta.SM,eta.SM)
  print(country.list[i])
  
}


#####-------------- Period 10 year --------------#####

db.SP <- db.MM <- db.SM <- 
  data.frame(vt=NA,eta=NA,eta.up=NA,eta.lo=NA,country=NA,
             rates=NA,rates.fit=NA,year=as.Date(NA),
             year.fit=as.Date(NA),model=NA)

for(i in 1:length(country.list)){ 
  db.SP <- rbind(db.SP,db.eu.SP[[i]]) 
  db.MM <- rbind(db.MM,db.eu.MM[[i]])
  db.SM <- rbind(db.SM,db.eu.SM[[i]])
}

db.SP <- db.SP[-1,]
db.MM <- db.MM[-1,]
db.SM <- db.SM[-1,]

db <- rbind(db.SP, db.MM, db.SM)


#####-------------- Save death file ------------#####

# save.image(file="fig/forecast/forecasts.EU.R.RData")
# load("fig/forecast/forecasts.EU.R.RData")

#####-------------- Themes fit and forecast --------------#####

type = c("Poisson Serfling"="SP", "Modulation model"="SP-STSS",
         "Smooth trend model"="SP-STFS")
col = c("darkgoldenrod3","cyan3","forestgreen")

fig.for <- function(data, country.name){
  
  data %>%
    filter(model %in% names(type), country == country.name) %>%
    mutate(type = factor(model, levels = names(type), labels = type)) %>%
    ggplot() +
    # Observed death counts
    geom_line(aes(x=year, y=log(rates)), color="red") +  
    geom_line(aes(x=year.fit, y=log(rates.fit)), color="grey") +
    # Predictions
    geom_line(aes(x=year, y=vt, colour=type), linetype="dashed") +
    geom_smooth(aes(x=year, y=eta, ymin=unname(eta.lo), ymax=unname(eta.up),
                    fill=type, colour=type), 
                size=0.5,linetype="solid", stat='identity') +
    facet_grid(vars(type), scales = "free") +
    scale_x_date(date_labels="%Y", date_breaks="2 year", expand=c(0.05,0.05)) +
    scale_fill_manual(values = col, labels = unname(type)) +
    scale_colour_manual(values = col, labels = unname(type)) +
    labs(title = country.name, x = "Month", y = "Log death rates") +
    theme_minimal() +
    theme(
      # Axis
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 10),
      # Facet
      strip.text.x = element_text(size = 10),
      strip.text.y = element_text(size = 10),
      # Legend
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.position="none", # "bottom"
      legend.key.width=unit(2.5, "cm"),
      # Title
      plot.title = element_text(hjust = 0.5, size = 10))
}

figEU63 <- fig.for(db, country.name="Sweden")
ggsave("fig/forecast/FigEU 63.png", plot=figEU63, width=6, height=7, bg="white")


######################################################
##########            AMPLITUDE             ##########
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

# Load function for forecasting
source("fun/ForFun total.R")

#####-------------- Modulation amplitude --------------#####

eu <- rates %>% filter(date>="2010-01-01", date<="2020-01-01")
period <- sort(unique(eu$date))
country.list=unique(eu$country)
# Start and end month of the forecasts in each window
tfit <- 121
# Weight vector (1 estimation set, 0 forecast set)
w <- rep(1,length(period))
# Number and degree of basis (10 years fit * 2 basis)
nx = 20
bdeg = 3

db <- as.data.frame(matrix(NA,length(1:tfit),7))
colnames(db) <- c("date","country","r","vt","Ddetrend","Ddeseaso","Ddeampli")
db$date <- as.Date(period[1:tfit])
ampl.country <- replicate(25, db, simplify=F)
names(ampl.country) <- country.list


for(i in 1:length(country.list)){
  
  d <- (eu %>% filter(country == country.list[i]) %>% arrange(date))$deaths
  e <- (eu %>% filter(country == country.list[i]) %>% arrange(date))$exposure
  x <- 1:length(d)
  d.fit <- d[w==1] # subset estimation dataset
  e.fit <- e[w==1]
  
  #####----------------- Grid search for lambdas -----------------######
  
  ##### Grid search for lambda1 and lambda2 (MM) and lambda (SM)
  lambda.seq <- 10^seq(4,7,length.out=30)
  mape.mat <- matrix(NA, length(lambda.seq), length(lambda.seq)) 
  rownames(mape.mat) <- colnames(mape.mat) <- lambda.seq
  mape.s.MM <- rep(NA,300)
  
  ##### Rolling window on the estimation set
  for (h in 1:length(lambda.seq)){
    for(g in 1:length(lambda.seq)){
      for(k in seq(49,109,12)){ 
        d.sel <- d.fit[1:(k+11)]
        e.sel <- e.fit[1:(k+11)]
        w.sel <- rep(1,k+11)
        w.sel[k:(k+11)] <- 0
        nx.sel <- (floor(k/12)+1)*2
        Output.MM <- PIRLS(d.sel, e.sel, w.sel, nx.sel, bdeg, 2, 1, 
                           lambda.seq[h], lambda.seq[g])$theta
        x.s <- 1:length(d.sel)
        X.s <- bspline(x.s, nx.sel, bdeg)
        eta.s.MM <- c(X.s$BB %*% Output.MM)[k:(k+11)]
        r.obs.s <- (d.fit/e.fit)[k:(k+11)]
        mape.s.MM[k] <- ((100/12)*sum(abs((r.obs.s-exp(eta.s.MM))/r.obs.s)))
      }
      mape.mat[h,g] <- mean(na.omit(mape.s.MM))
    }}
  mape.min.MM <- which(mape.mat == min(mape.mat), arr.ind = TRUE)
  lambda1 = lambda.seq[mape.min.MM[1]]
  lambda2 = lambda.seq[mape.min.MM[2]]
  
  #####----------------- Modulation model -----------------######
  
  theta.MM <- PIRLS(d, e, w, nx, bdeg, 2, 1, lambda1, lambda2)$theta
  X.MM <- bspline(x, nx, bdeg)$BB
  vt.MM <- c(X.MM[,1:(nx+bdeg)] %*% theta.MM[1:(nx+bdeg)])
  d <- d[1:tfit]
  e <- e[1:tfit]
  vt <- vt.MM[1:tfit]  
  
  # Compute the seasonal component and modulation amplitude
  nbas = nx+bdeg
  fCB <- X.MM[,(nbas+1):(2*(nbas))] %*% theta.MM[(nbas+1):(2*(nbas))]
  gSB <- X.MM[,(2*(nbas)+1):(3*(nbas))] %*% theta.MM[(2*(nbas)+1):(3*(nbas))]
  sumfg <- fCB + gSB
  fB <- X.MM[,1:(nbas)] %*% theta.MM[(nbas+1):(2*(nbas))]
  gB <- X.MM[,1:(nbas)] %*% theta.MM[(2*(nbas)+1):(3*(nbas))]
  rho <- sqrt(fB^2+gB^2)
  
  # Save the results
  ampl.country[[i]]$country <- country.list[i]
  ampl.country[[i]]$r <- d/e
  ampl.country[[i]]$vt <- vt
  ampl.country[[i]]$Ddetrend <- log((d/e)/exp(vt))
  ampl.country[[i]]$Ddeseaso <- sumfg[1:tfit]
  ampl.country[[i]]$Ddeampli <- rho[1:tfit]
  
  print(country.list[i])
  
}


#####-------------- Save death file ------------#####

# save.image(file="fig/forecast/forecasts.EU.R.modulation.RData")
# load("fig/forecast/forecasts.EU.R.amplitude.RData")


#####-------------- Plots fit and forecast --------------#####

ampl.db <- rbind(
  ampl.country[[1]], ampl.country[[2]], ampl.country[[3]],
  ampl.country[[4]], ampl.country[[5]], ampl.country[[6]],
  ampl.country[[7]], ampl.country[[8]], ampl.country[[9]],
  ampl.country[[10]], ampl.country[[11]], ampl.country[[12]],
  ampl.country[[13]], ampl.country[[14]], ampl.country[[15]],
  ampl.country[[16]], ampl.country[[17]], ampl.country[[18]],
  ampl.country[[19]], ampl.country[[20]], ampl.country[[21]],
  ampl.country[[22]], ampl.country[[23]], ampl.country[[24]],
  ampl.country[[25]])


figEU71 <- ampl.db %>% filter(country == "Sweden") %>%
  ggplot() +
  geom_line(aes(x=date, y=Ddetrend), color="grey", size=0.5) +  
  geom_line(aes(x=date, y=Ddeseaso), linetype="dashed") +
  geom_line(aes(x=date, y=Ddeampli), color="red", linetype="dashed", size=0.5) +
  geom_line(aes(x=date, y=-Ddeampli), color="red", linetype="dashed", size=0.5) +
  scale_x_date(date_labels="%Y", date_breaks="2 year", expand=c(0.05,0.05)) +  
  labs(title = "Sweden", x = "Month", y = "Log rates/trend") +
  theme_minimal() +
  theme(
    # Axis
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    # Legend
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position="none", # "bottom"
    legend.key.width=unit(2.5, "cm"),
    # Title
    plot.title = element_text(hjust = 0.5, size = 10))
ggsave("fig/forecast/FigEU 71.png", plot=figEU71, width=7, height=4, bg="white")

