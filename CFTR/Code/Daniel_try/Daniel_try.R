library(fgsea)
library(ggplot2)
library(scTenifoldNet)
library(Matrix)

load("CFTR/Daniel_results/SRS4245406.RData")
WT <- SRS4245406$WT
WT <- as.matrix(WT)
gKO = 'Cftr'

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
