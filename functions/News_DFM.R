# News_DFM()    Calculates changes in news
#
#  Syntax:
#  [y_old, y_new, singlenews, actual, fore, weight ,t_miss, v_miss, innov] = ...
#    News_DFM(X_old, X_new, Q, t_fcst, v_news) 
#
#  Description:
#   News DFM() inputs two datasets, DFM parameters, target time index, and 
#   target variable index. The function then produces Nowcast updates and
#   decomposes the changes into news.
#
#  Input Arguments:
#    X_old:  Old data matrix (old vintage)
#    X_new:  New data matrix (new vintage)
#    Res:    DFM() output results (see DFM for more details)
#    t_fcst: Index for target time
#    v_news: Index for target variable 
#
#  Output Arguments:
#    y_old:       Old nowcast
#    y_new:       New nowcast
#    single_news: News for each data series
#    actual:      Observed series release values
#    fore:        Forecasted series values
#    weight:      News weight
#    t_miss:      Time index for data releases
#    v_miss:      Series index for data releases
#    innov:       Difference between observed and predicted series values ("innovation")

News_DFM <- function(X_old, X_new, Res, t_fcst, v_news) {
  
  source("functions/para_const.R")
  library(purrr, quietly = T)
  library(pracma, quietly = T)
  
  # Initialize variables ----------------------------------------------------
  
  r <- ncol(Res$C)
  N <- ncol(X_new)
  singlenews <- matrix(0, 1, N) # Initialize news vector (will store news for each series)

  # No forecast case: Already values for variables v_news at time t_fcst --------

  if (!is.na(X_new[t_fcst, v_news])) {
    Res_old <- para_const(X_old, Res, 0) # Apply Kalman filter for old data
    
    for (i in 1:ncol(v_news)) { # Loop for each target variable
      
      # (Observed value) - (predicted value)
      singlenews[,v_news[i]] <- X_new[t_fcst, v_news[i]] - Res_old$X_sm[t_fcst, v_news[i]]
      
      # Set predicted and observed y values
      y_old[1,i] <- Res_old$X_sm[t_fcst, v_news[i]]
      y_new[1,i] <- X_new[t_fcst, v_news[i]]
      
    }
    
    # Forecast-related output set to empty
    actual <- NA
    forecast <- NA
    weight <- NA
    t_miss <- NA
    v_miss <- NA
    innov <- NA

  }
  
  # Forecast case (these are broken down into (A) and (B)) ------------------
  
  else {
    
    # Initialize series mean/standard deviation respectively
    Mx <- Res$Mx
    Wx = Res$Wx
    
    # Calculate indicators for missing values (1 if missing, 0 otherwise)
    miss_old <- is.na(X_old)
    miss_new <- is.na(X_new)
    
    # Indicator for missing--combine above information to single matrix where:
    # (i) -1: Value is in the old data, but missing in new data
    # (ii) 1: Value is in the new data, but missing in old data 
    # (iii) 0: Values are missing from/available in both datasets
    i_miss <- miss_old - miss_new
    
    # Time/variable indicies where case (b) is true
    t_miss <- which(i_miss == 1, arr.ind = T)[,1]
    v_miss <- which(i_miss == 1, arr.ind = T)[,2]
    
    # Forecast subcase (A): no new information --------------------------------

    if (is_empty(v_miss)) {
      
      # Fill in missing variables using a Kalman filter
      Res_old <- para_const(X_old, Res, 0)
      Res_new <- para_const(X_new, Res, 0)
      
      # Set predicted and observed y values. New y value is set to old
      y_old <- Res_old$X_sm[t_fcst, v_news]
      y_new <- y_old
      # y_new <- Res_new$X_sm(t_fcst, v_news)
      
      # No news, so nothing returned for news-related output
      groupnews <- NA
      singlenews <- NA
      gain <- NA
      gainSer <- NA
      actual <- NA
      forecast <- NA
      weight <- NA
      t_miss <- NA
      v_miss <- NA
      innov <- NA
      
    }
    
    # Forecast subcase (B): new information -----------------------------------
    
    else {

      # Difference between forecast time and new data time
      lag <- t_fcst - t_miss

      # Gives biggest time interval between forecast and new data
      k <- max(abs(lag), max(lag) - min(lag))
      
      C <- Res$C # Observation matrix
      R <- t(Res$R) # Covariance for observation matrix residuals
      
      # Number of new events
      n_news <- length(lag)
      
      # Smooth old dataset
      Res_old <- para_const(X_old, Res, k)
      Plag <- Res_old$Plag
            
      # Smooth new dataset
      Res_new <- para_const(X_new, Res, 0)
      
      # Subset for target variable and forecast time
      y_old <- Res_old$X_sm[t_fcst, v_news]
      y_new <- Res_new$X_sm[t_fcst, v_news]
      
      P <- Res_old$P[1:dim(Res_old$P)[1],1:dim(Res_old$P)[2],2:dim(Res_old$P)[3]]
      P1 <- NA # Initialize projection onto updates
      
      # Cycle through total number of updates
      for (i in 1:n_news) {
        
        h <- abs(t_fcst - t_miss[i])
        m <- max(t_miss[i], t_fcst)
        
        # If location of update is later than the forecasting date
        if (t_miss[i] > t_fcst) {
          Pp <- Plag[[h+1]][1:dim(Plag[[h+1]])[1],1:dim(Plag[[h+1]])[2],m] # P(1:r, h*r+1:h*r+r, m)'
        }
        else {
          Pp <- t(Plag[[h+1]][1:dim(Plag[[h+1]])[1],1:dim(Plag[[h+1]])[2],m]) # P(1:r, h*r+1:h*r+r, m)
        }
        P1 <- cbind(P1, Pp %*% C[v_miss[i], 1:r]) # Projection on updates
        
      }
      
      P1 <- P1[,-1]
      
      innov <- 1:length(t_miss)
      for (i in 1:length(t_miss)) {
        # Standardize predicted and observed values
        X_new_norm <- (X_new[t_miss[i], v_miss[i]] - Mx[v_miss[i]]) / Wx[v_miss[i]]
        X_sm_norm <- (Res_old$X_sm[t_miss[i], v_miss[i]] - Mx[v_miss[i]]) / Wx[v_miss[i]]
        
        # Innovation: Gives [observed] data - [predicted data]
        innov[i] <- X_new_norm - X_sm_norm
      }
      
      ins <- length(innov)
      P2 <- NA
      p2 <- NA
      WW <- matrix(0, nrow = max(v_miss), ncol = max(v_miss))
      
      # Gives non-standardized series weights
      for (i in 1:length(lag)) {
        
        for (j in 1:length(lag)) {
          
          h <- abs(lag[i] - lag[j])
          m <- max(t_miss[i], t_miss[j])
          
          if (t_miss[j] > t_miss[i]) {
            Pp <- Plag[[h+1]][1:dim(Plag[[h+1]])[1],1:dim(Plag[[h+1]])[2],m] # P(1:r,h*r+1:(h+1)*r,m)'
          }
          else {
            Pp <- t(Plag[[h+1]][1:dim(Plag[[h+1]])[1],1:dim(Plag[[h+1]])[2],m]) # P(1:r,h*r+1:(h+1)*r,m)
          }
          
          if ((v_miss[i] == v_miss[j]) && (t_miss[i] != t_miss[j])) {
            WW[v_miss[i], v_miss[j]] <- 0
          }
          else {
            WW[v_miss[i], v_miss[j]] <- R[v_miss[i], v_miss[j]]
          }
          
          p2 <- cbind(p2, C[v_miss[i], 1:r] %*% Pp %*% C[v_miss[j], 1:r] + WW[v_miss[i], v_miss[j]])
          
        }
        P2 <- rbind(P2, p2)
        p2 <- NA
      }
      P2 <- P2[-1,]
      P2 <- P2[,-1]
      
      
      totnews <- matrix(NA)
      temp <- array(NA, dim = c(1, n_news, 1))
      gain <- array(NA, dim = c(1, n_news, 1))
      for (i in 1:length(v_news)) { # loop on v_news
        # Convert to real units (unstadardized data)
        totnews[1, i] <- Wx[v_news[i]] %*% C[v_news[i], 1:r] %*% P1 %*% inv(P2) %*% innov
        temp[1:dim(temp)[1], 1:dim(temp)[2], i] <- Wx[v_news[i]] %*% C[v_news[i], 1:r] %*% P1 %*% inv(P2) * t(innov)
        gain[1:dim(gain)[1], 1:dim(gain)[2], i] <- Wx[v_news[i]] %*% C[v_news[i], 1:r] %*% P1 %*% inv(P2)
      }
      
      # Initialize output objects
      singlenews <- matrix(NA, max(t_miss) - min(t_miss) + 1, N) # ncol(v_news)
      actual <- matrix(NA, N, 1) # Actual forecasted values
      forecast <- matrix(NA, N, 1) # Forecasted values
      weight <- array(NA, dim = c(N, 1, length(v_news)))
      
      # Fill in output values
      for (i in 1:length(innov)) {
        actual[v_miss[i], 1] <- X_new[t_miss[i], v_miss[i]]
        forecast[v_miss[i], 1] <- Res_old$X_sm[t_miss[i], v_miss[i]]
        
        for (j in 1:length(v_news)) {
          singlenews[t_miss[i] - min(t_miss) + 1, v_miss[i]] <- temp[1, i, j]
          weight[v_miss[i], 1:dim(weight)[2], j] <- gain[1:dim(gain)[1], i, j] / Wx[v_miss[i]]
        }
      }
      
      singlenews <- colSums(singlenews) # Returns total news
      
      v_miss <- unique(v_miss)
      
    }
    
  }
  
  result <- list()
  result$y_old <- y_old
  result$y_new <- y_new
  result$singlenews <- singlenews
  result$actual <- actual
  result$forecast <- forecast
  result$weight <- weight
  result$t_miss <- t_miss
  result$v_miss <- v_miss
  result$innov <- innov
  result$ResOld <- Res_old
  result$ResNew <- Res_new
  
  return(result)
  
}
