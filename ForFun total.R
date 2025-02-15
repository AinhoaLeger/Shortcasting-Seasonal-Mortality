
#####-------------- Compute matrix B --------------#####

#     Input
#     X = x-domain vector
#     NDX = number of intervals in domain
#     BDEG = degree of B-spline (quadratic = 2, etc)
#     Output
#     B.MM = B-splines for trend and seasonality
#     B.SM = B-splines for trend

bspline <- function(X, NDX, BDEG){
  
  # Compute the proportion of the interval between every knot
  XL <- min(X)-1
  XR <- max(X)+1
  dx <- (XR - XL)/NDX

  # Place the knots
  knots <- seq(XL - BDEG*dx, XR + BDEG*dx, by=dx)
  
  # Evaluate the B-splines design matrix
  B <- spline.des(knots,X,BDEG+1,0*X)$design
  
  # Compute extended matrix for MM
  CB.MM <- diag(cos(((2*pi)/12)*X)) %*% B
  SB.MM <- diag(sin(((2*pi)/12)*X)) %*% B
  B.MM <- cbind(B,CB.MM,SB.MM)
  
  # Compute extended matrix for SM
  CB.SM <- cos(((2*pi)/12)*X)
  SB.SM <- sin(((2*pi)/12)*X)
  B.SM <- cbind(B,as.matrix(CB.SM),as.matrix(SB.SM)) 
  
  return(list(B=B, BB=B.MM, BB2=B.SM))
  
}

#####-------------- Compute penalty D'D --------------#####

#     Input
#     Same arguments as bspline().  Also
#     PORD = order of the penalty
#     LAMBDA = smoothing parameter
#     Output
#     P = penalty matrix

setup.P <- function(X, NDX, BDEG, PORD1, PORD2, LAMBDA1, LAMBDA2){

  # Evaluate the B-splines design matrix
  B <- bspline(X, NDX, BDEG)$B
  
  # Differences matrix for the trend
  D.T <- diff(diag(ncol(B)), diff=PORD1)
  DtD.T <- t(D.T) %*% D.T
  # Differences matrix for the seasonality
  D.S <- diff(diag(ncol(B)), diff=PORD2)
  DtD.S <- t(D.S) %*% D.S 
  
  # Penalty matrix - modulation model
  P.T <- kronecker(LAMBDA1,DtD.T)
  P.S <- kronecker(LAMBDA2,DtD.S)
  P.MM <- adiag(P.T,P.S,P.S)
  
  return(list(DPen=D.T, P=P.MM))
  
}

setup.P2 <- function(X, NDX, BDEG, PORD, LAMBDA){
  
  # Evaluate the B-splines design matrix
  B <- bspline(X, NDX, BDEG)$B
  
  # Differences matrix for the trend
  D.T <- diff(diag(ncol(B)), diff=PORD)
  DtD.T <- t(D.T) %*% D.T
  # Differences matrix for the fixed seasonality
  D.F <- diag(2)
  DtD.F <- t(D.F) %*% D.F 
  
  # Penalty matrix - smooth trend model
  P.T <- kronecker(LAMBDA,DtD.T)
  P.SM <- adiag(P.T,DtD.F)
  
  return(list(DPen=D.T, P2=P.SM))
  
}

#####-------------- Estimate the GLM Poisson --------------#####

#     Program to estimate regression parameters for a Poisson regression. 
#     Confidence intervals are included.
#     Main function: calculates
#          theta = estimated parameters
#          XtWX (needed for SE's)
#          Dev = deviance
#          Tr  = trace of hat matrix

IRLS <- function(X, V, D, E, THETA){
  
  # Iterations
  for(it in 1:20){
    Eta <- X %*% THETA
    Mu <- c(E*exp(Eta))
    W <- c(Mu*V)
    XtWX <- t(X) %*% (W*X)
    z <- Eta + (D-Mu)/Mu
    XtWz <- t(X) %*% (W*z)
    theta.old <- THETA
    THETA <- solve(XtWX, XtWz)
    tol <- max(abs(THETA - theta.old))/mean(abs(THETA))
    # cat(it, tol, "\n")
    if(tol < 10^(-8)) break
  }
  
  # Matrix variance-covariance
  VXt <- solve(XtWX)
  Var <- X %*% VXt %*% t(X)
  # Compute standard errors
  se.eta <- sqrt(diag(Var))
  
  # Penalized Poisson deviance
  y.init <- D
  y.init[D==0] <- 10^(-4) 
  # norm_vec <- function(x) sqrt(sum(x^2))
  Dev <- 2*sum( V*D*log(y.init/Mu) ) 
  # Hat matrix and trace
  H <- Var %*% W
  Tr <- sum(diag(H))
  # AIC and BIC criteria
  Aic <- Dev + 2*Tr
  Bic <- Dev + log( sum(V) )*Tr
  
  return(list(theta=THETA, se.eta=se.eta, Aic=Aic, Bic=Bic))
  
}

#####-------------- Estimate the modulation model --------------#####

#     Program to estimate regression parameters for a penalized 
#     Poisson regression. Confidence intervals are included.
#     Main function: calculates
#          tBWBpP = B'WB + P (needed for SE's)
#          Tr  = trace of hat matrix
#          Dev = deviance

PIRLS <- function(D, E, V, NDX, BDEG, PORD1, PORD2, LAMBDA1, LAMBDA2){
  
  # Extended design matrix
  x <- 1:length(D)
  B <- bspline(x, NDX, BDEG)$BB
  # Penalty
  P <- setup.P(x, NDX, BDEG, PORD1, PORD2, LAMBDA1, LAMBDA2)$P
  # Initialise theta
  DinitT <- matrix(log( D+1 ), nrow = nrow(B), ncol = 1)
  THETA <- solve( t(B) %*% B + P, t(B) %*% DinitT)
  # Iterations
  for(it in 1:20){
    Eta <- B %*% THETA
    Mu <- c(E*exp(Eta))
    W <- c(Mu*V)
    tBWB <- t(B) %*% (W*B)
    tBWBpP <- tBWB + P
    z <- Eta + (D-Mu)/Mu
    tBWz <- t(B)%*%(W*z)
    theta.old <- THETA
    THETA <- solve(tBWBpP, tBWz)
    tol <- max(abs(THETA - theta.old))/mean(abs(THETA))
    # cat(it, tol, "\n")
    if(tol < 10^(-8)) break
  }
  
  # Matrix variance-covariance
  VBt <- solve(tBWBpP) 
  Var <- B %*% VBt %*% t(B)  
  # Compute standard errors
  se.eta <- sqrt(diag(Var))
  
  # Penalized Poisson deviance
  y.init <- D
  y.init[D==0] <- 10^(-4) 
  D.pen <- setup.P(x, NDX, BDEG, PORD1, PORD2, LAMBDA1, LAMBDA2)$DPen
  alphas <- THETA[(1:(NDX+BDEG))]
  betas <- THETA[(1:(NDX+BDEG))+(NDX+BDEG)]
  gammas <- THETA[(1:(NDX+BDEG))+2*(NDX+BDEG)]
  Dev <- 2*sum( V* D*log(y.init/Mu) ) + 
    LAMBDA1 * t(alphas)%*%t(D.pen)%*%D.pen%*%alphas +
    LAMBDA2 * t(betas)%*%t(D.pen)%*%D.pen%*%betas + 
    LAMBDA2 * t(gammas)%*%t(D.pen)%*%D.pen%*%gammas
  # Hat matrix and trace
  H <- solve(tBWBpP, tBWB)
  Tr <- sum(diag(H))
  # AIC and BIC criteria
  Aic <- Dev + 2*Tr
  Bic <- Dev + log( sum(V) )*Tr
  
  return(list(theta=THETA, se.eta=se.eta, Aic=Aic, Bic=Bic))
  
}

#####-------------- Estimate the smooth trend model --------------#####

#     Program to estimate regression parameters for a penalized 
#     Poisson regression. Confidence intervals are included.
#     Main function: calculates
#          BtMBplusP = B'MB + P (needed for SE's)
#          Tr  = trace of hat matrix
#          Dev = deviance

PIRLS2 <- function(D, E, V, NDX, BDEG, PORD, LAMBDA){
  
  # Extended design matrix
  x <- 1:length(D)
  B <- bspline(x, NDX, BDEG)$BB2
  # Penalty
  P <- setup.P2(x, NDX, BDEG, PORD, LAMBDA)$P2
  # Initialise theta
  DinitT <- matrix(log( D+1 ), nrow = nrow(B), ncol = 1)
  THETA <- solve( t(B) %*% B + P, t(B) %*% DinitT)
  # Iterations
  for(it in 1:20){
    Eta <- B %*% THETA
    Mu <- c(E*exp(Eta))
    W <- c(Mu*V)
    tBWB <- t(B) %*% (W*B)
    tBWBpP <- tBWB + P
    z <- Eta + (D-Mu)/Mu
    tBWz <- t(B)%*%(W*z)
    theta.old <- THETA
    THETA <- solve(tBWBpP, tBWz)
    tol <- max(abs(THETA - theta.old))/mean(abs(THETA))
    # cat(it, tol, "\n")
    if(tol < 10^(-8)) break
  }
  
  # Matrix variance-covariance
  VBt <- solve(tBWBpP) 
  Var <- B %*% VBt %*% t(B)  
  # Compute standard errors
  se.eta <- sqrt(diag(Var))
  
  # Penalized Poisson deviance
  y.init <- D
  y.init[D==0] <- 10^(-4) 
  D.pen <- setup.P2(x, NDX, BDEG, PORD, LAMBDA)$DPen
  theta <- THETA[(1:(NDX+BDEG))]
  Dev <- 2*sum( V* D*log(y.init/Mu) ) + 
    LAMBDA * t(theta)%*%t(D.pen)%*%D.pen%*%theta
  # Hat matrix and trace
  H <- solve(tBWBpP, tBWB)
  Tr <- sum(diag(H))
  # AIC and BIC criteria
  Aic <- Dev + 2*Tr
  Bic <- Dev + log( sum(V) )*Tr
  
  return(list(theta=THETA, se.eta=se.eta, Aic=Aic, Bic=Bic))
  
}


