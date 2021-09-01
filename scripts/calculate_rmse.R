# Nowcast quality estimation ----------------------------------------------

load("output/rmse_perfomance/nwcst_res.RData")

nwcst_result %>% 
  dplyr::select(Date, nwcst_value, `Real GDP: Russia`) %>% 
  drop_na() %>% 
  mutate(nowcast = nwcst_value, fact = `Real GDP: Russia`) %>% 
  dplyr::select(Date, nowcast, fact) %>% 
  .[-1,] %>% 
  gather(key = "type", value = "value", nowcast:fact) %>% 
  ggplot(mapping = aes(x = Date, y = value, colour = type)) +
  geom_line()

nwcst_result %>% 
  dplyr::select(Date, nwcst_value, `Real GDP: Russia`) %>% 
  mutate(nowcast = nwcst_value, fact = `Real GDP: Russia`) %>% 
  dplyr::select(Date, nowcast, fact) %>% 
  fill(fact, .direction = "up") %>% 
  add_column(nwcst_type = rep(c("nowcast (1m)", "nowcast (2m)", "nowcast (3m)"), 20)) %>% 
  mutate(naive = dplyr::lag(fact, n = 3)) %>% 
  # mutate(naive_1m = dplyr::lag(fact, n = 4)) %>% 
  mutate(error_dfm = (fact - nowcast) ^ 2,
         error_naive = (fact - naive) ^ 2) %>% 
  .[-(1:3),] %>% 
  group_by(nwcst_type) %>% 
  summarise(relative_rmse = sum(error_dfm) / n() / (sum(error_naive) / n()))

# Forecast quality estimation ---------------------------------------------

load("output/rmse_perfomance/fcst_res.RData")

fcst_result %>% 
  dplyr::select(Date, fcst_value, `Real GDP: Russia`) %>% 
  drop_na() %>% 
  mutate(forecast = fcst_value, fact = `Real GDP: Russia`) %>% 
  dplyr::select(Date, forecast, fact) %>% 
  dplyr::filter(Date <= "2019-12-01") %>% 
  gather(key = "type", value = "value", forecast:fact) %>% 
  ggplot(mapping = aes(x = Date, y = value, colour = type)) +
  geom_line()

fcst_result %>% 
  dplyr::select(Date, fcst_value, `Real GDP: Russia`) %>% 
  mutate(forecast = fcst_value, fact = `Real GDP: Russia`) %>% 
  dplyr::select(Date, forecast, fact) %>% 
  fill(fact, .direction = "up") %>% 
  add_column(fcst_type = rep(c("forecast (1m)", "forecast (2m)", "forecast (3m)"), 20)) %>% 
  mutate(naive = dplyr::lag(fact, n = 6)) %>% 
  mutate(error_dfm = (fact - forecast) ^ 2,
         error_naive = (fact - naive) ^ 2) %>% 
  .[-(1:6),] %>% 
  group_by(fcst_type) %>% 
  summarise(relative_rmse = sum(error_dfm) / n() / (sum(error_naive) / n()))
