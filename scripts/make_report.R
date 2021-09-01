# Load functions ----------------------------------------------------------

source("functions/make_nowcasting_report.R")

# User inputs -------------------------------------------------------------

series <- "Real GDP: Russia" # Nowcasting "Real GDP: Russia"
period <- "2021Q3" # Nowcasting quarter
old <- "2021-08-30" # old vintage
new <- "2021-08-31" # new vintage
Res <- readRDS("output/ResDFM.rds") # Load estimated parameters

# Res$Mx[63] <- (1.015 ^ (0.25) - 1) * 100 # potential growth in MRC1

# Update nowcasting result ------------------------------------------------

make_nowcasting_report(series, period, old, new, Res)

