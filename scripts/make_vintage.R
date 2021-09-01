# Список с информацией о рядах --------------------------------------------

# Names of sheets to parse
sheet.names <- c("Real.modified", "Soft.modified", "World.modified", "Financial.modified", "Regional.modified")

spec_temp <- as_tibble() # датасет, в кот-м будет храниться информация о рядах
for (i in 1:length(sheet.names)) { # пробегаемся по всем листам и собираем данные
  data_temp <- "input/Data.File.xlsx" %>% 
    read_excel(sheet = sheet.names[i], col_names = F) %>% 
    .[(1:6),] %>% 
    t() %>% 
    as_tibble() %>% 
    .[-1, c(1, 3:6)]
  colnames(data_temp) <- c("SeriesRawName", "Source", "LastValue", "Updated", "Delay")
  data_temp <- data_temp %>% 
    mutate(Updated = as.Date(as.numeric(Updated), origin = "1899-12-30"),
           LastValue = as.Date(as.numeric(LastValue), origin = "1899-12-30"))
  spec_temp <- bind_rows(spec_temp, data_temp)
}

info_temp <- "input/Spec.File.xlsx" %>% 
  read_excel() %>% 
  dplyr::filter(Model == 1) %>% 
  dplyr::select(SeriesRawName, SeriesName, SA, Transformation)

spec_blocks <- "input/Spec.File.xlsx" %>% 
  read_excel() %>% 
  dplyr::filter(Model == 1) %>%
  .[,6:9] %>% 
  as.matrix()

spec_trans <- "input/Spec.File.xlsx" %>% 
  read_excel() %>% 
  dplyr::filter(Model == 1) %>%
  dplyr::select(Transformation) %>% 
  as.matrix()

spec_frequency <- "input/Spec.File.xlsx" %>% 
  read_excel() %>% 
  dplyr::filter(Model == 1) %>%
  dplyr::select(Frequency) %>% 
  as.matrix()

spec_base <- "input/Spec.File.xlsx" %>% 
  read_excel() %>% 
  dplyr::filter(Model == 1) %>% 
  dplyr::select(SeriesRawName, SeriesName, SeriesID, Category, Units) %>% 
  left_join(spec_temp, by = "SeriesRawName") %>% 
  mutate(Delay = as.numeric(Delay)) %>% 
  dplyr::select(SeriesName, SeriesID, Units, Category, Source, LastValue, Updated, Delay)

spec_delay <- spec_base %>% 
  dplyr::select(Delay) %>% 
  as.matrix()

spec_base <- spec_base %>% 
  dplyr::select(SeriesName, SeriesID, Units, Category, Source, LastValue, Updated)

info <- list(Base = spec_base, Blocks = spec_blocks, Frequency = spec_frequency,
             Transformation = spec_trans, Delay = spec_delay)
rm(spec_base, spec_blocks, spec_frequency, spec_trans, spec_temp, i, data_temp,
   spec_delay)

# Обработка данных --------------------------------------------------------

database <- as_tibble() # датасет, в кот-й будут записываться ряды
for (i in 1:length(sheet.names)) { # проходим по всем листам и собираем данные
  data_temp <- "input/Data.File.xlsx" %>% 
    read_excel(sheet = sheet.names[i]) %>% 
    .[-(1:5),] %>% 
    rename(Date = Indicator) %>% 
    mutate(Date = as.Date(as.numeric(Date), origin = "1899-12-30")) %>% 
    gather(key = "SeriesRawName", value = "Value", -Date) %>% 
    na_if(0) %>% 
    drop_na()
  database <- bind_rows(database, data_temp) # передаем данные в итоговый датасет
}

database <- database %>% 
  right_join(info_temp, by = "SeriesRawName") %>% 
  right_join(info$Base, by = "SeriesName") %>% 
  mutate(Value = as.numeric(Value)) %>% 
  dplyr::select(Date, SeriesName, Value, SA, Transformation, Category, Source, LastValue, Updated)

# Сглажка месячных индикаторов
data_for_sa <- database %>% 
  dplyr::filter(SA == 1) %>% 
  dplyr::filter(SeriesName != "Real GDP: Russia") # датасет рядов, кот-е нужно сгладить
series.names <- unique(data_for_sa$SeriesName)
data_new <- as_tibble()
for (i in 1:length(series.names)) {
  data_temp <- data_for_sa %>% 
    dplyr::filter(SeriesName == series.names[i])
  model <- seas(x = ts(data_temp$Value, start = c(year(min(data_temp$Date)), month(min(data_temp$Date))), frequency = 12))
  data_temp <- data_temp %>% 
    add_column(Value_SA = model$series$s11)
  data_new <- bind_rows(data_new, data_temp)
  print(series.names[i])
}

data_new <- data_new %>%
  mutate(SeriesName = as_factor(SeriesName)) %>%
  group_by(SeriesName) %>%
  mutate(Change_SA = if_else(Transformation == "pch", Value_SA / dplyr::lag(Value_SA) * 100 - 100,
                             if_else(Transformation == "chg", Value_SA - dplyr::lag(Value_SA),
                                     Value_SA)))

# Сглажка ВВП
data_temp <- database %>% 
  dplyr::filter(SeriesName == "Real GDP: Russia")
model <- seas(x = ts(data_temp$Value, start = c(year(min(data_temp$Date)), quarter(min(data_temp$Date))), frequency = 4),
              regression.variables = c("TC2020.2"), outlier = NULL)
# view(model)
data_temp <- data_temp %>% 
  add_column(Value_SA = c(model$series$s11)) %>% 
  mutate(Change_SA = if_else(Transformation == "pch", Value_SA / dplyr::lag(Value_SA) * 100 - 100,
                             if_else(Transformation == "chg", Value_SA - dplyr::lag(Value_SA),
                                     Value_SA)))
data_new <- bind_rows(data_new, data_temp)

# Соединяем с датасетом несглаживаемых рядов
Base <- database %>% 
  dplyr::filter(SA == 0) %>% 
  mutate(Value_SA = Value) %>% 
  mutate(SeriesName = as.factor(SeriesName)) %>%
  group_by(SeriesName) %>%
  mutate(Change_SA = if_else(Transformation == "pch", Value_SA / dplyr::lag(Value_SA) * 100 - 100,
                             if_else(Transformation == "chg", Value_SA - dplyr::lag(Value_SA),
                                     Value_SA))) %>%
  bind_rows(data_new) %>%
  dplyr::select(Date, SeriesName, Value_SA) %>% 
  spread(key = SeriesName, value = Value_SA)

X <- database %>% 
  dplyr::filter(SA == 0) %>% 
  mutate(Value_SA = Value) %>% 
  mutate(SeriesName = as.factor(SeriesName)) %>%
  group_by(SeriesName) %>%
  mutate(Change_SA = if_else(Transformation == "pch", Value_SA / dplyr::lag(Value_SA) * 100 - 100,
                             if_else(Transformation == "chg", Value_SA - dplyr::lag(Value_SA),
                                     Value_SA))) %>%
  bind_rows(data_new) %>%
  dplyr::select(Date, SeriesName, Change_SA) %>% 
  spread(key = SeriesName, value = Change_SA) %>%
  dplyr::select(-Date) %>%
  as.matrix() %>%
  ts(start = c(2002, 1), frequency = 12)

Time <- database %>%
  dplyr::select(Date) %>%
  as.matrix() %>% 
  as_date()

nwcst_dataset <- database %>% 
  dplyr::filter(SA == 0) %>% 
  mutate(Value_SA = Value) %>% 
  mutate(SeriesName = as.factor(SeriesName)) %>%
  group_by(SeriesName) %>%
  mutate(Change_SA = if_else(Transformation == "pch", Value_SA / dplyr::lag(Value_SA) * 100 - 100,
                             if_else(Transformation == "chg", Value_SA - dplyr::lag(Value_SA),
                                     Value_SA))) %>%
  bind_rows(data_new) %>%
  ungroup() %>% 
  dplyr::select(Date, SeriesName, Value, Value_SA, Change_SA, Category, Source, LastValue, Updated)

write.xlsx(x = nwcst_dataset, file = "output/data/temp_file.xlsx", col.names = T, row.names = T, sheetName = "data", append = F)

database <- list(Base = Base, X = X, Time = Time)

# Сохраняем результат -----------------------------------------------------

# Plot indicator
# data_df %>% 
#   dplyr::filter(SeriesName == "Real GDP: Russia") %>% 
#   dplyr::select(Date, Value_SA, SeriesName) %>% 
#   ggplot(mapping = aes(x = Date, y = Value_SA)) +
#   geom_line()

df <- list(database = database, info = info)

save(df, file = paste("input/vintages/", date(today()), ".RData", sep = ""))

# write.csv2(x = df$database$Base, file = "data.csv", sep = ";", dec = ".")
