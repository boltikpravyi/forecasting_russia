# SKF:Applies Kalman filter
# 
#  Syntax:
#    S <- SKF(Y, A, C, Q, R, Z_0, V_0)
#
#  Description:
#    SKF() applies the Kalman filter
#
#  Input parameters:
#    Y: k-by-nobs matrix of input data
#    A: m-by-m transition matrix 
#    C: k-by-m observation matrix
#    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
#    R: k-by-k covariance for observation matrix residuals (e_t)
#    Z_0: 1-by-m vector, initial value of state
#    V_0: m-by-m matrix, initial value of state covariance matrix
#
#  Output parameters:
#    S$Zm: m-by-nobs matrix, prior/predicted factor state vector
#          (S$Zm[,t] <- Z_t | t-1)
#    S$ZmU: m-by-(nobs+1) matrix, posterior/updated state vector
#           (S$Zm[t+1] <- Z_t | t)
#    S$Vm: m-by-m-by-nobs array, prior/predicted covariance of factor
#          state vector (S$Vm[,:,t] <- V_t | t-1)  
#    S$VmU: m-by-m-by-(nobs+1) array, posterior/updated covariance of
#           factor state vector (S$VmU[,,t+1] <- V_t | t)
#    S$loglik: scalar, value of likelihood function
#    S$k_t: k-by-m Kalman gain

SKF <- function(Y, A, C, Q, R, Z_0, V_0) {
  
  source("functions/MissData.R")
  library(purrr, quietly = T)
  
  # Initialize output values ------------------------------------------------
  
  # Output structure & dimensions of state space matrix
  m <- ncol(C)
  
  # Outputs time for data matrix. "number of observations"
  nobs <- ncol(Y)
  
  # Instantiate output
  S <- list()
  S$Zm <- matrix(NA, m, nobs) # Z_t | t-1 (prior)
  S$Vm <- array(NA, dim = c(m, m, nobs)) # V_t | t-1 (prior)
  S$ZmU <- matrix(NA, m, (nobs+1)) # Z_t | t (posterior/updated)
  S$VmU <- array(NA, dim = c(m, m, (nobs+1))) # V_t | t (posterior/updated)
  S$loglik <- 0
  
  # Set initial values ------------------------------------------------------
  
  Zu <- Z_0 # Z_0 | 0 (In below loop, Zu gives Z_t | t)
  Vu <- V_0 # V_0 | 0 (In below loop, Vu guvse V_t | t)
  
  # Store initial values
  S$ZmU[,1] <- Zu
  S$VmU[,,1] <- Vu
  
  # Kalman filter procedure -------------------------------------------------
  
  for (t in 1:nobs) {
    
    ### Calculating prior distribution
    
    # Use transition eqn to create prior estimate for factor
    # i.e. Z = Z_t | t-1
    Z <- A %*% Zu
    
    # Prior covariance matrix of Z (i.e. V = V_t | t-1)
    #   Var(Z) = Var(A*Z + u_t) = Var(A*Z) + Var(\epsilon) = 
    #   A*Vu*A' + Q
    
    V <- (A %*% Vu %*% t(A)) + Q
    V <- 0.5 * (V + t(V)) # Trick to make symmetric
    
    ### Calculating posterior distribution
    
    # Removes missing series: These are removed from Y, C, and R
    data <- MissData(Y[,t], C, R)
    Y_t <- data$y
    C_t <- data$C
    R_t <- data$R
    
    # Check if y_t contains no data. If so, replace Zu and Vu with prior.
    if (is_empty(Y_t)) {
      Zu <- Z
      Vu <- V
    }
    else {
      # Steps for variance and population regression coefficients:
      # Var(c_t * Z_t + e_t) = c_t Var(A) c_t' + Var(u) = c_t*V *c_t' + R
      VC <- V %*% t(C_t)  
      iF <- solve((C_t %*% VC) + R_t)
      
      # Matrix of population regression coefficients (QuantEcon eqn #4)
      VCF <- VC %*% iF
      
      # Gives difference between actual and predicted observation matrix values
      innov <- Y_t - (C_t %*% Z)
      
      # Update estimate of factor values (posterior)
      Zu <- Z  + (VCF %*% innov)
      
      # Update covariance matrix (posterior) for time t
      Vu <- V  - (VCF %*% t(VC))
      Vu <- 0.5 * (Vu + t(Vu)) # Approximation trick to make symmetric
      
      # Update log likelihood 
      S$loglik <- S$loglik + 0.5 * (log(det(iF)) - t(innov) %*% iF %*% innov)
    }
    
    ### Store output
    
    # Store covariance and observation values for t-1 (priors)
    S$Zm[,t] <- Z
    S$Vm[,,t] <- V
    
    # Store covariance and state values for t (posteriors)
    # i.e. Zu = Z_t | t & Vu = V_t | t
    S$ZmU[,(t+1)] <- Zu
    S$VmU[,,(t+1)] <- Vu
    
  }
  
  # Store Kalman gain k_t
  if (is_empty(Y_t)) {
    S$k_t <- matrix(0, m, m)
  }
  else {
    S$k_t <- VCF %*% C_t
  }
  
  return(S)
  
}
