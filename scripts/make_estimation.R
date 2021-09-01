###################### Dynamic factor model (DFM) #########################
########### This script estimates a dynamic factor model (DFM) ############
############# using a panel of monthly and quarterly series. ##############

# Load functions ----------------------------------------------------------
source("functions/dfm.R")

# Load model specification and dataset ------------------------------------
load("input/vintages/2021-08-31.RData")
df$database$X <- df$database$X[,c(1:49, 51:63, 50)]
vintage_thrshld <- "2014-12-01"
vintage <- PRTDB(mts = df$database$X, delay = df$info$Delay, vintage = vintage_thrshld)
X <- ts(vintage[-(1:3),], start = c(2002, 4), frequency = 12)
X <- ts(X[-nrow(X),], start = c(2002, 4), frequency = 12)
Spec <- list(SeriesID = as.matrix(df$info$Base[,1]), SeriesName = as.matrix(df$info$Base[,1]),
             Frequency = df$info$Frequency, Blocks = df$info$Blocks, BlockNames = colnames(df$info$Blocks))

# Run dynamic factor model (DFM) and save estimation output ---------------
threshold <- 1e-5 # Set to 1e-5 for more robust estimates
Res <- dfm(X, Spec, threshold)
saveRDS(Res, file = "Res.rds")

# Display output ----------------------------------------------------------

# Table with monthly factor loadings
{nQ <- sum(str_count(string = Spec$Frequency, pattern = "q")) # Number of quarterly series
nM <- sum(str_count(string = Spec$Frequency, pattern = "m")) # Number of monthly series
nLags <- 5 # 5 comes from monthly-quarterly aggregation
nFactors <- sum(matrix(1, 1, ncol(Spec$Blocks)))

monthly_loadings <- as.tibble(Res$C[(1:nM), seq(1, nFactors*5, 5)])
colnames(monthly_loadings) <- Spec$BlockNames
monthly_loadings <- monthly_loadings %>% 
  add_column(Indicator = Spec$SeriesName[1:nM], .before = "Global")
options(pillar.sigfig = 3)
print(monthly_loadings, n = nrow(monthly_loadings))}

# Plot GDP and its Common Factor ------------------------------------------

idx <- 30
GlobalFactor <- Res$Z[,1:5] %*% Res$C[idx, 1:5]
SoftFactor <- Res$Z[,6:11] %*% Res$C[idx, 6:11]
RealFactor <- Res$Z[,11:16] %*% Res$C[idx, 11:16]
FinancialFactor <- Res$Z[,16:21] %*% Res$C[idx, 16:21]
RegionalFactor <- Res$Z[,21:26] %*% Res$C[idx, 21:26]

CommonComponent <- GlobalFactor + RealFactor + RegionalFactor + SoftFactor
plot(Res$x_sm[,idx], type = "l")
lines(CommonComponent, type = "p")

# Plot monthly series and its Common Factor -------------------------------

idx <- 12
GlobalFactor <- Res$Z[,1] * Res$C[idx, 1]
SoftFactor <- Res$Z[,6] * Res$C[idx, 6]
RealFactor <- Res$Z[,11] * Res$C[idx, 11]
FinancialFactor <- Res$Z[,16] * Res$C[idx, 16]
RegionalFactor <- Res$Z[,21] * Res$C[idx, 21]

CommonComponent <- GlobalFactor + RealFactor + RegionalFactor + SoftFactor
plot(Res$X_sm[,idx], type = "l")
lines(CommonComponent, type = "p")
