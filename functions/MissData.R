# Syntax:
# Description:
#   Eliminates the rows in y & matrices C, R that correspond to missing 
#   data (NA) in y
#
# Input:
#   y: Vector of observations at time t
#   C: Observation matrix
#   R: Covariance for observation matrix residuals
#
# Output:
#   y: Vector of observations at time t (reduced)     
#   C: Observation matrix (reduced)     
#   R: Covariance for observation matrix residuals
#   L: Used to restore standard dimensions(n x #) where # is the nr of 
#      available data in y

MissData <- function(y, C, R) {
  
  # Returns 1 for nonmissing series
  ix <- !is.na(y)
  
  # Index for columns with nonmissing variables
  e <- diag(length(y))
  L <- e[,ix]
  
  # Removes missing series
  y <- y[ix]
  
  # Removes missing series from observation matrix
  if (sum(ix) == 1) {
    C <- t(C[ix,])
  }
  else {
    C <- C[ix,]
  }

  # Removes missing series from transition matrix
  R <- R[ix,ix]
  
  return(list(y = y, C = C, R = R, L = L))
  
}
