update_nowcast <- function(X_old, X_new, Time, Spec, Res, series, period, vintage_old, vintage_new) {
    
    source("functions/News_DFM.R")
    source("functions/para_const.R")
    source("functions/runKF.R")
    source("functions/MissData.R")
    
    vintage_old <- lubridate::as_date(vintage_old)
    vintage_new <- lubridate::as_date(vintage_new)
    
    # Make sure datasets are the same size
    N <- ncol(X_new)
    T_old <- nrow(X_old)
    T_new <- nrow(X_new)
    if (T_new > T_old) {
        X_old <- rbind(X_old, matrix(NA, T_new - T_old, N))
    }
    # append 1 year (12 months) of data to each dataset to allow for forecasting at different horizons
    X_old <- rbind(X_old, matrix(NA, 3, N))
    X_new <- rbind(X_new, matrix(NA, 3, N))
    
    y <- lubridate::year(Time[length(Time)])
    m <- lubridate::month(Time[length(Time)])
    d <- lubridate::day(Time[length(Time)])
    
    Time <- as.Date(union(Time, (as.Date(str_c(y, "-", m, "-", d)) + months(1:3))))
    
    i_series <- which(series == Spec$SeriesID)
    
    series_name <- Spec$SeriesName[i_series]
    freq <- Spec$Frequency[i_series]
    
    switch(freq,
           "M" = {
               y <- period
               m <- freq
               d <- 1
               t_nowcast = which(Time == as.Date(y, m, d))
           },
           "Q" = {
               y <- as.numeric(str_extract(period, "\\d{4}"))
               q <- as.numeric(str_extract(str_extract(period, "Q\\d{1}"), "\\d{1}"))
               m <- 3 * q
               d <- 1
               t_nowcast <- which(Time == as.Date(str_c(y, "-", m, "-", d)))
           }
    )
    
    if (is_empty(t_nowcast)) {
        stop("Period is out of nowcasting horizon (up to one year ahead)")
    }
    
    # Update nowcast for target variable 'series' (i) at horizon 'period' (t)
    #   > Relate nowcast update into news from data releases:
    #     a. Compute the impact from data revisions
    #     b. Compute the impact from new data releases
    
    X_rev <- X_new
    X_rev[is.na(X_old)] <- NA

    # Compute news ------------------------------------------------------------
    
    # Compute impact from data revisions
    result_revision <- News_DFM(X_old, X_rev, Res, t_nowcast, i_series)

    # Compute impact from data releases
    result_release <- News_DFM(X_rev, X_new, Res, t_nowcast, i_series)
    
    # Summarize change in model forecast
    cat('Table: Model forecast decomposition \n')
    
    if (sum(!is.na(result_release$actual)) != 0) {
        nowcast_impacts <- tryCatch(
            {
                tibble(series_name = Spec$SeriesID[result_release$v_miss]) %>%
                    mutate(forecast = na.omit(result_release$forecast[1:length(result_release$forecast)]),
                           actual = na.omit(result_release$actual[1:length(result_release$actual)])) %>%
                    mutate(news = actual - forecast) %>%
                    mutate(weight = na.omit(result_release$weight[1:length(result_release$weight)])) %>%
                    mutate(impact = weight * news)
            }
        )
        print(nowcast_impacts, n = nrow(nowcast_impacts))
    }
    
    results <- list()
    results$release <- result_release
    results$revision <- result_revision
    results$nowcast_impacts <- nowcast_impacts
    return(results)
}
