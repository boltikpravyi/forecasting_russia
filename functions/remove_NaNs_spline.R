# Description:
# 
# remNaNs    Treats NaNs in dataset for use in DFM.
# 
#   Syntax:
#     [X,indNaN] = remNaNs(X,options)
# 
#   Description:
#     remNaNs() processes NaNs in a data matrix X according to 5 cases (see
#     below for details). These are useful for running functions in the 
#     'DFM.m' file that do not take missing value inputs.
# 
#   Input parameters:
#     X (T x n): Input data where T gives time and n gives the series. 
#     options: A structure with two elements:
#       options.method (numeric):
#       - 1: Replaces all missing values using filter().
#       - 2: Replaces missing values after removing trailing and leading
#            zeros (a row is 'missing' if >80% is NaN)
#       - 3: Only removes rows with leading and closing zeros
#       - 4: Replaces missing values after removing trailing and leading
#            zeros (a row is 'missing' if all are NaN)
#       - 5: Replaces missing values with spline() then runs filter().
# 
#       options.k (numeric): remNaNs() relies on MATLAB's filter function
#       for the 1-D filter. k controls the rational transfer function
#       argument's numerator (the denominator is set to 1). More
#       specifically, the numerator takes the form 'ones(2*k+1,1)/(2*k+1)'
#       For additional help, see MATLAB's documentation for filter().
# 
#   Output parameters:
#     X: Outputted data. 
#     indNaN: A matrix indicating the location for missing values (1 for NaN).

remNaNs_spline <- function(X, options) {
  
  library(pracma, quietly = T)
  
  t <- nrow(X) # Gives dimensions for data input
  N <- ncol(X) # Gives dimensions for data input
  k <- options$k # Inputted options
  indNaN <- is.na(X) # Returns location of NAs
  
  switch(options$method,
         "1" = { # replace all the missing values
           for (i in 1:N) { # loop through columns
             x <- X[,i]
             x[indNaN[,i]] <- median(x, na.rm = T) # Replace missing values series median 
             x_MA <- as.matrix(signal::filter(matrix(1, 2*k + 1, 1) / (2*k + 1), 1, rbind(as.matrix(x[1] * matrix(1, k, 1)), as.matrix(x), as.matrix(x[length(x)] * matrix(1, k, 1)))))
             x_MA <- x_MA[(2*k + 1):length(x_MA)] # Match dimensions
             # Replace missing observations with filtered values
             x[indNaN[,i]] <- x_MA[indNaN[,i]]
             X[,i] <- x # Replace vector
           }
         },
         "2" = { # replace missing values after removing leading and closing zeros
           # Returns row sum for NaN values. Marks true for rows with more than 80% NaN
           rem1 <- rowSums(indNaN) > N * 0.8
           nanLead <- cumsum(rem1) == seq(1, t)
           nanEnd <- cumsum(rev(rem1)) == seq(1, t)
           nanEnd <- rev(nanEnd) # Reverses nanEnd
           nanLE <- nanLead | nanEnd
           
           # Subsets X for for 
           X <- X[!nanLE,]
           indNaN <- is.na(X) # Index for missing values
           # Loop for each series
           for (i in 1:N) {
             x <- X[,i]
             isnanx <- is.na(x)
             t1 <- min(which(!isnanx)) # First non-NaN entry 
             t2 <- max(which(!isnanx)) # Last non-NaN entry
             # Interpolates without NaN entries in beginning and end
             x[t1:t2] <- interp1(x = which(!isnanx), y = x[!isnanx], xi = seq(t1, t2), method = "spline")
             isnanx <- is.na(x)
             # replace NaN observations with median
             x[isnanx] <- median(x, na.rm = T)
             # Apply filter
             x_MA <- as.matrix(signal::filter(matrix(1, 2*k + 1, 1) / (2*k + 1), 1, rbind(as.matrix(x[1] * matrix(1, k, 1)), as.matrix(x), as.matrix(x[length(x)] * matrix(1, k, 1)))))
             x_MA <- x_MA[(2*k + 1):length(x_MA)] # Match dimensions
             # Replace nanx with filtered observations
             x[isnanx] <- x_MA[isnanx]
             X[,i] <- x
           }
         },
         "3" = { # only remove rows with leading and closing zeros
           rem1 <- (rowSums(indNaN) == N)
           nanLead <- (cumsum(rem1) == seq(1, t))
           nanEnd <- (cumsum(rev(rem1)) == seq(1, t))
           nanEnd <- rev(nanEnd) # Reverses nanEnd
           nanLE <- nanLead | nanEnd
           X <- X[!nanLE,]
           indNaN <- is.na(X) # Index for missing values
         },
         "4" = { # remove rows with leading and closing zeros & replace missing values
           rem1 <- rowSums(indNaN) == N
           nanLead <- cumsum(rem1) == seq(1, t)
           nanEnd <- cumsum(rev(rem1)) == seq(1, t)
           nanEnd <- rev(nanEnd) # Reverses nanEnd
           nanLE <- nanLead | nanEnd
           X <- X[-nanLE,]
           indNaN <- is.na(X) # Index for missing values
           for (i in 1:N) {
             x <- X[,i]
             isnanx <- is.na(x)
             t1 <- min(which(!isnanx)) # First non-NaN entry
             t2 <- max(which(!isnanx)) # Last non-NaN entry
             # Interpolates without NaN entries in beginning and end
             x[t1:t2] <- interp1(x = which(!isnanx), y = x[!isnanx], xi = seq(t1, t2), method = "spline")
             isnanx <- is.na(x)
             # replace NaN observations with median
             x[isnanx] <- median(x, na.rm = T)
             # Apply filter
             x_MA <- as.matrix(signal::filter(matrix(1, 2*k + 1, 1) / (2*k + 1), 1, rbind(as.matrix(x[1] * matrix(1, k, 1)), as.matrix(x), as.matrix(x[length(x)] * matrix(1, k, 1)))))
             x_MA <- x_MA[(2*k + 1):length(x_MA)] # Match dimensions
             # Replace nanx with filtered observations
             x[isnanx] <- x_MA[isnanx]
             X[,i] <- x
           }
         },
         "5" = { # replace missing values
           indNaN <- is.na(X)
           for (i in 1:N) {
             x <- X[,i]
             isnanx <- is.na(x)
             t1 <- min(which(!isnanx)) # First non-NaN entry
             t2 <- max(which(!isnanx)) # Last non-NaN entry
             # Interpolates without NaN entries in beginning and end
             x[t1:t2] <- spline(x = which(!isnanx), y = x[!isnanx], seq(t1, t2))$y
             isnanx <- is.na(x)
             # replace NaN observations with median
             x[isnanx] <- median(x, na.rm = T)
             # Apply filter
             x_MA <- as.matrix(signal::filter(matrix(1, 2*k + 1, 1) / (2*k + 1), 1, rbind(as.matrix(x[1] * matrix(1, k, 1)), as.matrix(x), as.matrix(x[length(x)] * matrix(1, k, 1)))))
             x_MA <- x_MA[(2*k + 1):length(x_MA)] # Match dimensions
             # Replace nanx with filtered observations
             x[isnanx] <- x_MA[isnanx]
             X[,i] <- x
           }
         },
         {
           print("nothing")
         }
  )
  return(list(X = X, indNaN = indNaN))
}
