# para_const()    Implements Kalman filter for "News_DFM.m"
#
#  Syntax:
#    Res <- para_const(X, P, lag)
#
#  Description:
#    para_const() implements the Kalman filter for the news calculation
#    step. This procedure smooths and fills in missing data for a given 
#    data matrix X. In contrast to runKF(), this function is used when
#    model parameters are already estimated.
#
#  Input parameters:
#    X: Data matrix. 
#    P: Parameters from the dynamic factor model.
#    lag: Number of lags
#
#  Output parameters:
#    Res [struc]: A structure containing the following:
#      Res$Plag: Smoothed factor covariance for transition matrix
#      Res$P:    Smoothed factor covariance matrix
#      Res$X_sm: Smoothed data matrix
#      Res$F:    Smoothed factors
#
# Kalman filter with specified paramaters written for 
# "MAXIMUM LIKELIHOOD ESTIMATION OF FACTOR MODELS ON DATA SETS WITH ARBITRARY
# PATTERN OF MISSING DATA." by Marta Banbura and Michele Modugno

para_const <- function(X, P, lag) {

  source("functions/SKF.R")
  source("functions/FIS.R")
  library(pracma, quietly = T)

  # Set model parameters and data preparation -------------------------------

  # Set model parameters
  Z_0 <- P$Z_0
  V_0 <- P$V_0
  A <- P$A
  C <- P$C
  Q <- P$Q
  R <- P$R
  Mx <- P$Mx
  Wx <- P$Wx

  # Prepare data
  TT <- nrow(X)

  # Standardise x
  Y <- t((X - kronecker(matrix(1, TT, 1), t(Mx))) / kronecker(matrix(1, TT, 1), t(Wx)))

  # Apply Kalman filter and smoother ----------------------------------------
  # See runKF() for details about FIS and SKF
  
  Sf <- SKF(Y, A, C, Q, R, Z_0, V_0) # Kalman filter
  
  Ss <- FIS(A, Sf) # Smoothing step
  
  # Calculate parameter output ----------------------------------------------

  Vs <- Ss$VmT[,,2:dim(Ss$VmT)[3]] # Smoothed factor covariance for transition matrix
  Vf <- Sf$VmU[,,2:dim(Sf$VmU)[3]] # Filtered factor posterior covariance
  Zsmooth <- Ss$ZmT # Smoothed factors
  Vsmooth <- Ss$VmT # Smoothed covariance values
  
  Plag <- list()
  Plag[[1]] <- Vs
  
  if (lag != 0) {
    for (jk in 1:lag) {
      Plag[[jk+1]] <- array(0, dim = c(nrow(Vf), ncol(Vf), dim(Plag[[1]])[3]))
      for (jt in dim(Plag[[1]])[3]:(lag + 1)) {
        As <- Vf[1:dim(Vf)[1],1:dim(Vf)[2],(jt-jk)] %*% t(A) %*% pinv(A %*% Vf[1:dim(Vf)[1],1:dim(Vf)[2],(jt - jk)] %*% t(A) + Q)
        Plag[[jk+1]][1:dim(Plag[[jk+1]])[1],1:dim(Plag[[jk+1]])[2],jt] <- As %*% Plag[[jk]][1:dim(Plag[[jk]])[1],1:dim(Plag[[jk]])[2],jt]
      }
    }
  }
  
  # Prepare data for output
  Zsmooth <- t(Zsmooth)
  x_sm <- Zsmooth[2:nrow(Zsmooth),] %*% t(C) # Factors to series representation
  X_sm <- (kronecker(matrix(1, TT, 1), t(Wx)) * x_sm) + kronecker(matrix(1, TT, 1), t(Mx)) # Standardized to unstandardized

  # Loading the structure with the results
  Res <- list()
  Res$Plag <- Plag
  Res$P <- Vsmooth
  Res$X_sm <- X_sm  
  Res$FF <- Zsmooth[2:nrow(Zsmooth),]
  
  return(Res)
  
}


  

