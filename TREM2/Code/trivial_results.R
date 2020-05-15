## check with trivial method
load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[which(g_wei > 0.2)]
source("TREM2/Code/utility.R")
E <- enrichFunction(gList, fdrThreshold = 0.05, nCategories = 20)
length(unique(E$Term))
