# Load functions ----------------------------------------------------------

source("functions/update_nowcast.R")

# User inputs -------------------------------------------------------------

series <- "Real GDP: Russia" # Nowcasting Real GDP: Russia
period <- "2021Q3" # Nowcasting quarter
Res <- readRDS("ResDFM.rds") # Load estimated parameters
# Res$Mx[63] <- (1.015 ^ (0.25) - 1) * 100 # potential growth in MRC1

# Construct old vintage
load("input/vintages/2021-08-31.RData") # load old dataset
df$database$X <- df$database$X[,c(1:49, 51:63, 50)]
vintage_old <- "2021-08-27"
vintage <- PRTDB(mts = df$database$X, delay = df$info$Delay, vintage = vintage_old)
X_old <- ts(vintage[-(1:3),], start = c(2002, 4), frequency = 12)
X_old <- ts(X_old[-nrow(X_old),], start = c(2002, 4), frequency = 12)

# Construct new vintage
load("input/vintages/2021-08-31.RData") # load new dataset
df$database$X <- df$database$X[,c(1:49, 51:63, 50)]
vintage_new <- "2021-08-31"
vintage <- PRTDB(mts = df$database$X, delay = df$info$Delay, vintage = vintage_new)
X_new <- ts(vintage[-(1:3),], start = c(2002, 4), frequency = 12)
X_new <- ts(X_new[-nrow(X_new),], start = c(2002, 4), frequency = 12)

# Update nowcast and decompose changes into news --------------------------

Time <- df$database$Base %>% 
  dplyr::filter(Date < vintage_new) %>% 
  dplyr::select(Date) %>% 
  .[-(1:3),] %>% 
  as.matrix() %>% 
  as_date()

Spec <- list(SeriesID = as.matrix(df$info$Base[,1]), SeriesName = as.matrix(df$info$Base[,1]),
             Frequency = df$info$Frequency)

results <- update_nowcast(X_old, X_new, Time, Spec, Res, series, period, vintage_old, vintage_new)

nowcast_impacts <- results$nowcast_impacts %>% 
  rename(SeriesName = series_name) %>% 
  left_join(df$info$Base, by = "SeriesName") %>% 
  dplyr::select(SeriesName, forecast, actual, news, weight, impact, Category, Source, LastValue, Updated) %>% 
  rename(ReleaseDate = Updated) %>% 
  add_column(Updated = rep(vintage_new, nrow(results$nowcast_impacts)))

nowcast_value <- tibble(Growth = results$release$y_new, Updated = vintage_new)

if (sum(!str_detect(dir("output/"), pattern = paste0("nowcasting_data_", period, ".RData"))) == 0) {
  nowcasting_data_2021Q3 <- list(impacts = nowcast_impacts, value = nowcast_value)
  save(nowcasting_data_2021Q3, file = "output/nowcasting_data_2021Q3.RData")
} else {
  load(paste0("output/nowcasting_data_", period, ".RData"))
  nowcasting_data_2021Q3$impacts <- bind_rows(nowcasting_data_2021Q3$impacts, nowcast_impacts)
  nowcasting_data_2021Q3$value <- bind_rows(nowcasting_data_2021Q3$value, nowcast_value)
  save(nowcasting_data_2021Q3, file = "output/nowcasting_data_2021Q3.RData")
}

write.csv2(x = nowcasting_data_2021Q3$value, file = "nowcasting_data.csv", sep = ";", dec = ".")
write.csv2(x = nowcasting_data_2021Q3$impacts, file = "nowcasting_data_imp.csv", sep = ";", dec = ".")
