library(fgsea)
library(ggplot2)
library(scTenifoldNet)
library(Matrix)

load("HNF4A-HNF4G/Daniel_Results/GSM3477499.RData")
WT <- GSM3477499$WT
WT <- as.matrix(WT)
gKO = c('Hnf4a','Hnf4g')

## row as 0
KO <- WT
KO[gKO,] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
saveRDS(MA, "Daniel_ori_row0.rds")

## column as 0
KO <- WT
KO[, gKO] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
saveRDS(MA, "Daniel_ori_col0.rds")

## both as 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
saveRDS(MA, "Daniel_ori_both0.rds")
