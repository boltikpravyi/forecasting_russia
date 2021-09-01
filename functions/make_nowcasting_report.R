make_nowcasting_report <- function(series, period, old, new, Res) {
  
  source("functions/update_nowcast.R") # load function for nowcast updating 
  
  # Old vintage
  load(paste0("input/vintages/", old, ".RData")) # load old dataset
  X_old <- ts(df$database$X[-(1:3),c(1:49, 51:63, 50)], start = c(2002, 4), frequency = 12)
  
  # New vintage
  load(paste0("input/vintages/", new, ".RData")) # load new dataset
  X_new <- ts(df$database$X[-(1:3),c(1:49, 51:63, 50)], start = c(2002, 4), frequency = 12)
  
  # Update nowcast and decompose changes into news --------------------------
  
  Time <- df$database$Base %>% 
    dplyr::select(Date) %>% 
    .[-(1:3),] %>% 
    as.matrix() %>% 
    as_date()
  
  Spec <- list(SeriesID = as.matrix(df$info$Base[,1]), SeriesName = as.matrix(df$info$Base[,1]),
               Frequency = df$info$Frequency)
  
  results <- update_nowcast(X_old, X_new, Time, Spec, Res, series, period, old, new)
  
  nowcast_impacts <- results$nowcast_impacts %>% 
    rename(SeriesName = series_name) %>% 
    left_join(df$info$Base, by = "SeriesName") %>% 
    dplyr::select(SeriesName, Units, forecast, actual, news, weight, impact, Category, Source, LastValue, Updated) %>% 
    rename(ReleaseDate = Updated) %>% 
    add_column(Updated = rep(new, nrow(results$nowcast_impacts))) %>% 
    mutate(Updated = as_date(Updated))
  
  nowcast_value <- tibble(Growth = results$release$y_new, Updated = new) %>% 
    mutate(Updated = as_date(Updated))
  
  write.xlsx(x = nowcast_impacts, file = "output/report/temp_file.xlsx", col.names = T,
             row.names = T, sheetName = "Impacts", append = F)
  write.xlsx(x = nowcast_value, file = "output/report/temp_file.xlsx", col.names = T,
             row.names = T, sheetName = "Value", append = T)
  
  
  
}