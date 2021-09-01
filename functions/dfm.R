# dfm: Runs the dynamic factor model
#
# Syntax:
#   Res <- DFM(X, Par)
#
# Description:
#   DFM() inputs the organized and transformed data X and parameter structure Par.xThen, the function outputs dynamic factor model
#   structure Res and data summary statistics (mean and standard deviation).
#
# Input arguments:
#   X: Kalman-smoothed data where missing values are replaced by their expectation
#   Par: A structure containing the following parameters:
#     $ blocks: Block loadings.
#     $ nQ: Number of quarterly series
#     $ p: Number of lags in transition matrix
#     $ r: Number of common factors for each block
#
# Output Arguments:
#
#   Res - structure of model results with the following fields
#     $ X_sm | Kalman-smoothed data where missing values are replaced by their expectation
#     $ Z | Smoothed states. Rows give time, and columns are organized according to Res.C.
#     $ C | Observation matrix. The rows correspond to each series, and the columns are organized as shown below:
#        - 1-20: These columns give the factor loadings. For example, 1-5 give loadings for the first block and are organized in
#        reverse-chronological order (f^G_t, f^G_t-1, f^G_t-2, f^G_t-3, f^G_t-4). Columns 6-10, 11-15, and 16-20 give loadings for
#        the second, third, and fourth blocks respectively.
#     $ R: Covariance for observation matrix residuals
#     $ A: Transition matrix. This is a square matrix that follows the same organization scheme as Res.C's columns. Identity matrices
#         are used to account for matching terms on the left and righthand side. For example, we place an I4 matrix to account for
#         matching (f_t-1; f_t-2; f_t-3; f_t-4) terms.
#     $ Q: Covariance for transition equation residuals.
#     $ Mx: Series mean
#     $ Wx: Series standard deviation
#     $ Z_0: Initial value of state
#     $ V_0: Initial value of covariance matrix
#     $ r: Number of common factors for each block
#     $ p: Number of lags in transition equation
#
# References:
#
# Marta Banbura, Domenico Giannone and Lucrezia Reichlin
# Nowcasting (2010)
# Michael P. Clements and David F. Hendry, editors,
# Oxford Handbook on Economic Forecasting.

dfm <- function(X, Spec, threshold) {
  
  source("functions/InitCond.R")
  source("functions/remove_NaNs_spline.R")
  source("functions/EMstep.R")
  source("functions/runKF.R")
  source("functions/em_converged.R")
  
  # Store model parameters --------------------------------------------------
  
  # DFM input specifications: See documentation for details
  Par <- list()
  Par$blocks <- Spec$Blocks # Block loading structure
  Par$nQ <- sum(str_count(string = Spec$Frequency, pattern = "q")) # Number of quarterly series
  Par$p <- 1 # Number of lags in autoregressive of factor (same for all factors)
  Par$r <- matrix(1, 1, ncol(Spec$Blocks)) # Number of common factors for each block
  
  table <- as.tibble(cbind(Spec$SeriesName, Spec$Blocks))
  colnames(table) <- cbind("Indicator", t(Spec$BlockNames))
  
  cat("Table 3: Block Loading Structure\n")
  print(table, n = nrow(table))
  
  cat("Estimating the dynamic factor model (DFM) ... \n\n\n") # Print a message to the console
  
  TT <- nrow(X)
  N <- ncol(X)
  r <- Par$r
  p <- Par$p
  nQ <- Par$nQ
  blocks <- Par$blocks
  dims <- dim(blocks)
  blocks <- as.numeric(blocks)
  dim(blocks) <- dims
  
  i_idio <- rbind(matrix(1, (N - nQ), 1), matrix(0, nQ, 1)) > 0
  
  # R * Lambda = q # Constraints on the loading of the quarterly variables
  
  R_mat <- rbind(c(2, -1, 0, 0, 0),
                 c(3, 0, -1, 0, 0),
                 c(2, 0, 0, -1, 0),
                 c(1, 0, 0, 0, -1))
  
  q <- matrix(0, 4, 1)
  
  if (missing(threshold)) {
    threshold <- 1e-5 # EM loop threshold (default value)
  }
  
  max_iter <- 1000 # EM loop maximum number of iterations
  
  # Prepare data ------------------------------------------------------------
  
  Mx <- apply(X, 2, mean, na.rm = T)
  Wx <- apply(X, 2, sd, na.rm = T)
  
  xNaN <- (X - kronecker(matrix(1, TT, 1), t(Mx))) / kronecker(matrix(1, TT, 1), t(Wx)) # Standardize series
  
  # Initial Conditions ------------------------------------------------------
  
  optNaN <- list()
  optNaN$method <- 2 # Remove leading and closing zeros
  optNaN$k <- 3 # Setting for filter(): See remNaN_spline
  
  data <- InitCond(xNaN, r, p, blocks, optNaN, R_mat, q, nQ, i_idio)

  A <- data$A
  C <- data$C
  Q <- data$Q
  R <- data$R
  Z_0 <- data$Z_0
  V_0 <- data$V_0
  
  # Initialize EM loop values
  previous_loglik <- -Inf
  num_iter <- 0
  LL <- -Inf
  converged <- 0
  
  # y for the estimation is WITH missing data
  y <- t(xNaN)

  # EM Loop -----------------------------------------------------------------
  
  # The model can be written as
  # y = C*Z + e;
  # Z = A*Z(-1) + v
  # where y is NxT, Z is (pr)xT, etc
  
  # Remove the leading and ending nans
  optNaN$method <- 3
  y_est <- t(remNaNs_spline(X = xNaN, options = optNaN)$X)
  
  while ((num_iter < max_iter) & !converged) { # Loop until converges or max iter
    # Applying EM algorithm
    data <- EMstep(y_est, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ, i_idio, blocks)
    C_new <- data$C_new
    R_new <- data$R_new
    A_new <- data$A_new
    Q_new <- data$Q_new
    Z_0 <- data$Z_0
    V_0 <- data$V_0
    loglik <- data$loglik
    
    C <- C_new
    R <- R_new
    A <- A_new
    Q <- Q_new
    
    if (num_iter > 2) { # Checking convergence
      data <- em_converged(loglik, previous_loglik, threshold, 1)
      converged <- data$converged
      decrease <- data$decrease[num_iter+1]
    }
    
    if ((num_iter %% 5 == 0) && (num_iter > 0)) { # Print updates to command window
      cat(str_c("Now running the ", num_iter, "th iteration of max ", max_iter, "\n"))
      cat("  LogLik", "         (%) Change", "\n")
      cat("  ", sprintf("%.3f", loglik), "        (", sprintf("%.3f", (100 * (loglik - previous_loglik) / previous_loglik)), "%)", "\n\n", sep = "")
    }
    
    LL <- list(LL = LL, loglik = loglik)
    previous_loglik <- loglik
    num_iter <- num_iter + 1
  }
  
  if (num_iter < max_iter) {
    cat("Successful: Convergence at ", num_iter, " iterations")
    
  }
  else {
    cat("Stopped because maximum iterations reached")
  }
  
  # Final run of the Kalman filter
  Zsmooth <- t(runKF(y, A, C, Q, R, Z_0, V_0)$zsmooth)

  x_sm <- Zsmooth[2:nrow(Zsmooth),] %*% t(C) # Get smoothed X
  
  # Loading the structure with the results ----------------------------------

  Res <- list()
  Res$x_sm <- x_sm
  
  Res$X_sm <- (kronecker(matrix(1, TT, 1), t(Wx)) * x_sm) + kronecker(matrix(1, TT, 1), t(Mx)) # Unstandardized, smoothed
  Res$Z <- Zsmooth[2:nrow(Zsmooth),]
  Res$C <- C
  Res$R <- R
  Res$A <- A
  Res$Q <- Q
  Res$Mx <- Mx
  Res$Wx <- Wx
  Res$Z_0 <- Z_0
  Res$V_0 <- V_0
  Res$r <- r
  Res$p <- p
  
  return(Res)
  
}
