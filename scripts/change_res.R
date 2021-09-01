ResDFM <- readMat("archive/ResDFM.mat")

names(ResDFM$Res) <- c("x_sm", "X_sm", "Z", "C", "R", "A", "Q", "Mx", "Wx", "Z_0", "V_0", "r", "p")
ResDFM <- ResDFM$Res
ResDFM <- as.list(ResDFM[,,1])
names(ResDFM) <- c("x_sm", "X_sm", "Z", "C", "R", "A", "Q", "Mx", "Wx", "Z_0", "V_0", "r", "p")
ResDFM$Mx <- c(ResDFM$Mx)
ResDFM$Wx <- c(ResDFM$Wx)
ResDFM$Z_0 <- c(ResDFM$Z_0)
ResDFM$p <- c(ResDFM$p)

saveRDS(ResDFM, file = "output/ResDFM.rds")
