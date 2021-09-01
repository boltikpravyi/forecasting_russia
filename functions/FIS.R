# FIS:    Applies fixed-interval smoother
#
#  Syntax:
#    S <- FIS(A, S)
#
#  Description:
#    SKF() applies a fixed-interval smoother, and is used in conjunction 
#    with SKF(). See  page 154 of 'Forecasting, structural time series models 
#    and the Kalman filter' for more details (Harvey, 1990).
#
#  Input parameters:
#    A: m-by-m transition matrix 
#    S: structure returned by SKF()
#
#  Output parameters:
#    S: FIS() adds the following smoothed estimates to the S structure: 
#    - S$ZmT: m-by-(nobs+1) matrix, smoothed states
#             (S$ZmT[,t+1] = Z_t|T) 
#    - S$VmT: m-by-m-by-(nobs+1) array, smoothed factor covariance
#             matrices (S$VmT[,,t+1] = V_t|T = Cov(Z_t|T))
#    - S$VmT_1: m-by-m-by-nobs array, smoothed lag 1 factor covariance
#               matrices (S$VmT_1[,,t] = Cov(Z_t Z_t-1|T))
#
#  Model:
#   Y_t = C_t Z_t + e_t for e_t ~ N(0, R)
#   Z_t = A Z_{t-1} + mu_t for mu_t ~ N(0, Q)

FIS <- function(A, S) {
  
  library(pracma, quietly = T)
  library(MASS, quietly = T)
  
  # Organize imput ----------------------------------------------------------
  
  # Initialize output matrices
  m <- nrow(S$Zm)
  nobs <- ncol(S$Zm)  
  S$ZmT <- matrix(0, m, (nobs+1))
  S$VmT <- array(0, dim = c(m, m, (nobs+1)))
  
  # Fill the final period of ZmT, VmT with SKF() posterior values
  S$ZmT[,(nobs+1)] <- drop(S$ZmU[,(nobs+1)])
  S$VmT[,,(nobs+1)] <- drop(S$VmU[,,(nobs+1)])
  
  # Initialize VmT_1 lag 1 covariance matrix for final period
  S$VmT_1 <- array(0, dim = c(m, m, nobs))
  S$VmT_1[,,nobs] <- (diag(m) - S$k_t) %*% A %*% drop(S$VmU[,,nobs])
  
  # Used for recursion process. See companion file for details
  J_2 <- drop(S$VmU[,,nobs]) %*% t(A) %*% pinv(drop(S$Vm[,,nobs]))
  
  # Run smoothing algorithm -------------------------------------------------
  
  # Loop through time reverse-chronologically (starting at final period nobs)
  
  for (t in nobs:1) {
    
    # Store posterior and prior factor covariance values
    VmU <- drop(S$VmU[,,t])
    Vm1 <- drop(S$Vm[,,t])
    
    # Store previous period smoothed factor covariance and lag-1 covariance
    V_T <- drop(S$VmT[,,(t+1)])
    V_T1 <- drop(S$VmT_1[,,t])
    
    J_1 <- J_2
    
    # Update smoothed factor estimate
    S$ZmT[,t] <- S$ZmU[,t] + (J_1 %*% (S$ZmT[,(t+1)] - (A %*% S$ZmU[,t])))
    
    # Update smoothed factor covariance matrix
    S$VmT[,,t] <- VmU + (J_1 %*% (V_T - Vm1) %*% t(J_1))
    
    if (t > 1) {
      # Update weight
      J_2 <- drop(S$VmU[,,(t-1)]) %*% t(A) %*% pinv(drop(S$Vm[,,(t-1)]), tol = .Machine$double.eps^(3))
      
      # Update lag 1 factor covariance matrix
      S$VmT_1[,,(t-1)] <- (VmU %*% t(J_2)) + (J_1 %*% (V_T1 - (A %*% VmU))) %*% t(J_2)
    }
    
  }
  
  return(S)
  
}
