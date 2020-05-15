## check with trivial method
load('NKX2-1/Daniel_Results/GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)
gKO <- 'Nkx2-1'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[which(g_wei > 0.4)]
source("NKX2-1/Code/utility.R")
E <- enrichFunction(gList, fdrThreshold = 0.05, nCategories = 20)
length(unique(E$Term))

