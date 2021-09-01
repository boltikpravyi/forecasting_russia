source("functions/update_nowcast.R")

Res <- readRDS("output/ResDFM.rds") # Load estimated parameters

# Set a threshold (for estimation)
inpt <- tibble(vintage_old = seq(ymd("2015-01-28"), ymd("2019-12-28"), by = "1 month"),
               vintage_new = seq(ymd("2015-02-28"), ymd("2020-01-28"), by = "1 month"),
               fcst_period = c(rep("2015q2", 3), rep("2015q3", 3), rep("2015q4", 3), rep("2016q1", 3),
                               rep("2016q2", 3), rep("2016q3", 3), rep("2016q4", 3), rep("2017q1", 3),
                               rep("2017q2", 3), rep("2017q3", 3), rep("2017q4", 3), rep("2018q1", 3),
                               rep("2018q2", 3), rep("2018q3", 3), rep("2018q4", 3), rep("2019q1", 3),
                               rep("2019q2", 3), rep("2019q3", 3), rep("2019q4", 3), rep("2020q1", 3)))

series <- "Real GDP: Russia" # Forecasting Real GDP: Russia

# User inputs
load("input/vintages/2021-08-31.RData")
df$database$X <- df$database$X[,c(1:49, 51:63, 50)]

Spec <- list(SeriesID = as.matrix(df$info$Base[,1]), SeriesName = as.matrix(df$info$Base[,1]),
             Frequency = df$info$Frequency)

fcst_result <- as_tibble()

for (i in 1:nrow(inpt)) {
  
  period <- inpt$fcst_period[i] # Nowcasting quarter
  
  # Construct old vintage
  vintage <- PRTDB(mts = df$database$X, delay = df$info$Delay, vintage = inpt$vintage_old[i])
  X_old <- ts(vintage[-(1:3),], start = c(2002, 4), frequency = 12)
  X_old <- ts(X_old[-nrow(X_old),], start = c(2002, 4), frequency = 12)
  
  # Construct new vintage
  vintage <- PRTDB(mts = df$database$X, delay = df$info$Delay, vintage = inpt$vintage_new[i])
  X_new <- ts(vintage[-(1:3),], start = c(2002, 4), frequency = 12)
  X_new <- ts(X_new[-nrow(X_new),], start = c(2002, 4), frequency = 12)
  
  Time <- df$database$Base %>% 
    dplyr::filter(Date < inpt$vintage_new[i]) %>% 
    dplyr::select(Date) %>% 
    .[-(1:3),] %>% 
    as.matrix() %>% 
    as_date()
  
  results <- update_nowcast(X_old, X_new, Time, Spec, Res, series, period, inpt$vintage_old[i], inpt$vintage_new[i])
  
  temp_result <- tibble(fcst_period = inpt$fcst_period[i], fcst_value = results$release$y_new)
  
  fcst_result <- bind_rows(fcst_result, temp_result)
  
  print(temp_result)
  
}

gdp_fact <- df$database$Base %>% 
  dplyr::select(Date, `Real GDP: Russia`) %>% 
  dplyr::filter(Date >= "2015-04-01") %>% 
  dplyr::filter(Date <= "2020-03-01")

fcst_result <- bind_cols(fcst_result, gdp_fact)

save(x = fcst_result, file = "output/rmse_perfomance/fcst_res.RData")
