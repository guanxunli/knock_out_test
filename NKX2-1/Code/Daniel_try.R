library(fgsea)
library(ggplot2)
library(scTenifoldNet)
library(Matrix)


load("../Results/GSM3716703.RData")
WT <- GSM3716703$WT
WT <- as.matrix(WT)
gKO <- 'Nkx2-1'

## row as 0
KO <- WT
KO[gKO,] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
set.seed(1)
DR <- scTenifoldNet::dRegulation(MA)
outputList <- list()
outputList$WT <- WT
outputList$KO <- KO
outputList$MA <- MA
outputList$diffRegulation <- DR
saveRDS(outputList, "Daniel_ori_row0.rds")

## column as 0
KO <- WT
KO[, gKO] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
set.seed(1)
DR <- scTenifoldNet::dRegulation(MA)
outputList <- list()
outputList$WT <- WT
outputList$KO <- KO
outputList$MA <- MA
outputList$diffRegulation <- DR
saveRDS(outputList, "Daniel_ori_col0.rds")

## both as 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
set.seed(1)
DR <- scTenifoldNet::dRegulation(MA)
outputList <- list()
outputList$WT <- WT
outputList$KO <- KO
outputList$MA <- MA
outputList$diffRegulation <- DR
saveRDS(outputList, "Daniel_ori_both0.rds")

##################################### sccnet try ############################################
WT <- readRDS("network_results/net_try.rds")
WT <- as.matrix(WT)
gKO <- 'Nkx2-1'

## row as 0
KO <- WT
KO[gKO,] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
set.seed(1)
DR <- scTenifoldNet::dRegulation(MA)
outputList <- list()
outputList$WT <- WT
outputList$KO <- KO
outputList$MA <- MA
outputList$diffRegulation <- DR
saveRDS(outputList, "Daniel_scc_row0.rds")

## column as 0
KO <- WT
KO[, gKO] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
set.seed(1)
DR <- scTenifoldNet::dRegulation(MA)
outputList <- list()
outputList$WT <- WT
outputList$KO <- KO
outputList$MA <- MA
outputList$diffRegulation <- DR
saveRDS(outputList, "Daniel_scc_col0.rds")

## both as 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(WT, KO)
set.seed(1)
DR <- scTenifoldNet::dRegulation(MA)
outputList <- list()
outputList$WT <- WT
outputList$KO <- KO
outputList$MA <- MA
outputList$diffRegulation <- DR
saveRDS(outputList, "Daniel_scc_both0.rds")
