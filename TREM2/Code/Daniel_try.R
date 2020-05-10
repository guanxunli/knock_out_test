library(fgsea)
library(ggplot2)
library(scTenifoldNet)
library(Matrix)

load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'

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
