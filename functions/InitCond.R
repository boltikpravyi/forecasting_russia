# InitCond: Calculates initial conditions for parameter estimation
#
#  Description:
#    Given standardized data and model information, InitCond() creates
#    initial parameter estimates. These are intial inputs in the EM
#    algorithm, which re-estimates these parameters using Kalman filtering
#    techniques.
#
# Inputs:
#  - x:      Standardized data
#  - r:      Number of common factors for each block
#  - p:      Number of lags in transition equation
#  - blocks: Gives series loadings
#  - optNaN: Option for missing values in spline. See remNaNs_spline() for details.
#  - Rcon:   Incorporates estimation for quarterly series (i.e. "tent structure")
#  - q:      Constraints on loadings for quarterly variables
#  - NQ:     Number of quarterly variables
#  - i_idio: Logical. Gives index for monthly variables (1) and quarterly (0)
#
# Output:
#  - A:   Transition matrix
#  - C:   Observation matrix
#  - Q:   Covariance for transition equation residuals
#  - R:   Covariance for observation equation residuals
#  - Z_0: Initial value of state
#  - V_0: Initial value of covariance matrix

InitCond <- function(x, r, p, blocks, optNaN, Rcon, q, nQ, i_idio) {
  
  library(pracma, quietly = T)
  source("functions/remove_NaNs_spline.R")
  
  pC <- ncol(Rcon) # Gives "tent" structure size (quarterly to monthly)
  ppC <- max(p, pC)
  n_b <- ncol(blocks) # Number of blocks
  
  data <- remNaNs_spline(x, optNaN) # Spline without NaNs
  xBal <- data$X
  indNaN <- data$indNaN
  
  TT <- nrow(xBal) # Time T series number N
  N <- ncol(xBal) # Time T series number N
  nM <- N - nQ # Number of monthly series
  
  xNaN <- xBal
  xNaN[indNaN] <- NA # Set missing values equal to NA
  res <- xBal # Spline output equal to res: Later this is used for residuals
  resNaN <- xNaN # Later used for residuals
  
  # Initialize model coefficient output
  C <- matrix(NA)
  A <- matrix(NA)
  Q <- matrix(NA)
  V_0 <- matrix(NA)
  
  # Set the first observations as NAs: For quarterly-monthly aggreg. scheme
  indNaN[1:(pC-1),] <- T
  
  for (i in 1:n_b) { # Loop for each block
    
    r_i <- r[i] # r_i = 1 when block is loaded
    
    ### Observation equation
    
    C_i <- matrix(0, N, r_i * ppC) # Initialize state variable matrix helper
    idx_i <- which(blocks[,i] == 1) # Returns series index loading block i
    idx_iM <- idx_i[idx_i < (nM + 1)] # Monthly series indicies for loaded blocks
    idx_iQ <- idx_i[idx_i > nM] # Quarterly series indicies for loaded blocks
    
    # Returns eigenvector v w/largest eigenvalue d
    data <- eigen(cov(res[,idx_iM]))
    v <- data$vectors[,1]
    d <- data$values[1]
    
    # Flip sign for cleaner output. Gives equivalent results without this section
    if (sum(v) < 0) {
      v <- -v
    }
    
    # For monthly series with loaded blocks (rows), replace with eigenvector
    # This gives the loading
    C_i[idx_iM, 1:r_i] <- v
    f <- res[,idx_iM] %*% v # Data projection for eigenvector direction
    FF <- matrix(NA)
    
    # Lag matrix using loading. This is later used for quarterly series
    for (kk in 1:((max(p + 1, pC) - 1) + 1)) {
      if (kk == 1) {
        FF <- f[(pC-kk+1):(length(f)-kk+1),]
      }
      else {
        FF <- cbind(FF, f[(pC-kk+1):(length(f)-kk+1),])
      }
    }
    
    Rcon_i <- kronecker(Rcon, diag(r_i)) # Quarterly-monthly aggregation scheme
    q_i <- kronecker(q, matrix(0, r_i, 1))
    
    # Produces projected data with lag structure (so pC - 1 fewer entries)
    ff <- FF[, 1:(r_i*pC)]
    
    for (j in t(idx_iQ)) { # Loop for quarterly variables
      
      # For series j, values are dropped to accommodate lag structure
      xx_j <- resNaN[pC:nrow(resNaN), j]
      
      if (sum(!is.na(xx_j)) < (ncol(ff) + 2)) {
        xx_j <- res[pC:nrow(res), j] # Replaces xx_j with spline if too many NaNs
      }
      
      ff_j <- ff[!is.na(xx_j),]
      xx_j <- xx_j[!is.na(xx_j)]
      
      iff_j <- solve(t(ff_j) %*% ff_j)
      Cc <- iff_j %*% t(ff_j) %*% xx_j # Least squares
      
      # Spline data monthly to quarterly conversion
      Cc <- Cc - (iff_j %*% t(Rcon_i) %*% solve(Rcon_i %*% iff_j %*% t(Rcon_i)) %*% ((Rcon_i %*% Cc) - q_i))
      
      C_i[j, 1:(pC*r_i)] <- t(Cc) # Place in output matrix
      
    }
    
    ff <- rbind(matrix(0, pC-1, pC*r_i), ff) # Zeros in first pC-1 entries (replace dropped from lag)
    
    # Residual calculations
    res <- res - (ff %*% t(C_i))
    resNaN <- res
    resNaN[indNaN] <- NA
    
    if (is.na(C)) {
      C <- C_i # Combine past loadings together
    }
    else {
      C <- cbind(C, C_i)
    }
    
    ### Transition equation
    
    z <- FF[,(1:r_i)] # Projected data (no lag)
    Z <- FF[,(r_i+1):(r_i*(p+1))] # Data with lag 1
    
    A_i <- t(matrix(0, r_i*ppC, r_i*ppC)) # Initialize transition matrix
    
    A_temp <- solve(t(Z) %*% Z) %*% t(Z) %*% z # OLS: gives coefficient value AR(p) process
    A_i[(1:r_i), 1:(r_i*p)] <- t(A_temp)
    A_i[(r_i+1):nrow(A_i), 1:(r_i*(ppC-1))] <- diag(r_i*(ppC-1))
    
    ####
    Q_i <- matrix(0, ppC*r_i, ppC*r_i)
    e <- z - (Z %*% A_temp) # VAR residuals
    Q_i[1:r_i, 1:r_i] <- cov(e) # VAR covariance matrix
    
    initV_i <- matrix((solve(diag((r_i*ppC)^2) - kronecker(A_i, A_i)) %*% c(Q_i)), r_i*ppC, r_i*ppC)
    
    # Gives top left block for the transition matrix
    if (is.na(A) && is.na(Q) && is.na(V_0)) {
      A <- A_i
      Q <- Q_i
      V_0 <- initV_i
    }
    else {
      A <- blkdiag(A, A_i)
      Q <- blkdiag(Q, Q_i)
      V_0 <- blkdiag(V_0, initV_i)
    }
    
  }
  
  eyeN <- diag(N) # Used inside observation matrix
  eyeN <- eyeN[,t(i_idio)]
  
  C <- cbind(C, eyeN)
  C <- cbind(C, rbind(matrix(0, nM, 5*nQ), t(kronecker(diag(nQ), c(1, 2, 3, 2, 1)))))
  R <- diag(apply(resNaN, 2, var, na.rm = T)) # Initialize covariance matrix for transition matrix
  
  ii_idio <- which(i_idio == 1) # Indicies for monthly variables
  n_idio <- length(ii_idio) # Number of monthly variables
  BM <- matrix(0, n_idio, n_idio) # Initialize monthly transition matrix values
  SM <- matrix(0, n_idio, n_idio) # Initialize monthly residual covariance matrix values
  
  for (i in 1:n_idio) { # Loop for monthly variables
    
    # Set observation equation residual covariance matrix diagonal
    R[ii_idio[i], ii_idio[i]] <- 1e-04
    
    # Subsetting series residuals for series i
    res_i <- resNaN[,ii_idio[i]]
    
    # Returns number of leading/ending zeros
    leadZero <- max(which(t(1:TT) == cumsum(is.na(res_i))))
    if (max(which(t(1:TT) == cumsum(is.na(rev(res_i))))) == -Inf) {
      endZero <- 0
    }
    else {
      endZero <- max(which(t(1:TT) == cumsum(is.na(rev(res_i)))))
    }
    
    # Truncate leading and ending zeros
    res_i <- res[,ii_idio[i]]
    res_i <- res_i[-((length(res_i) - endZero) + (1:length(res_i)))]
    res_i <- res_i[-(1:leadZero)]
    
    # Linear regression: AR 1 process for monthly series residuals
    BM[i,i] <- solve(t(res_i[1:(length(res_i)-1)]) %*% res_i[1:(length(res_i)-1)]) %*% t(res_i[1:(length(res_i)-1)]) %*% res_i[2:length(res_i)]
    SM[i,i] <- var(res_i[2:length(res_i)] - (res_i[1:(length(res_i)-1)] * BM[i,i])) # Residual covariance matrix
    
  }
  
  Rdiag <- diag(R)
  sig_e <- Rdiag[(nM+1):N] / 19
  Rdiag[(nM+1):N] <- 1e-04
  R <- diag(Rdiag) # Covariance for obs matrix residuals
  
  # For BQ, SQ
  rho0 <- 0.1
  temp <- matrix(0, 5, 5)
  temp[1,1] <- 1
  
  # Blocks for covariance matrices
  
  if (length(sig_e) == 1) {
    SQ <- kronecker((1-(rho0^2))*sig_e, temp)
    BQ <- kronecker(diag(nQ), rbind(cbind(rho0, matrix(0, 1, 4)), cbind(diag(4), matrix(0, 4, 1))))
    initViQ <- matrix(solve(diag((5*nQ)^2) - kronecker(BQ, BQ)) %*% c(SQ), 5*nQ, 5*nQ)
    initViM <- diag(1 / diag(diag(nrow(BM)) - BM^2)) %*% SM
  }
  else {
    SQ <- kronecker(diag((1-(rho0^2))*sig_e), temp)
    BQ <- kronecker(diag(nQ), rbind(cbind(rho0, matrix(0, 1, 4)), cbind(diag(4), matrix(0, 4, 1))))
    initViQ <- matrix(solve(diag((5*nQ)^2) - kronecker(BQ, BQ)) %*% c(SQ), 5*nQ, 5*nQ)
    initViM <- diag(1 / diag(diag(nrow(BM)) - BM^2)) %*% SM
  }
  
  # Output
  A <- blkdiag(A, BM, BQ) # Observation matrix
  Q <- blkdiag(Q, SM, SQ) # Residual covariance matrix (transition)
  Z_0 <- matrix(0, nrow(A), 1) # States
  V_0 <- blkdiag(V_0, initViM, initViQ) # Covariance of states
  
  return(list(A = A, C = C, Q = Q, R = R, Z_0 = Z_0, V_0 = V_0))
  
}
