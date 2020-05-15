## check with trivial method
load("CFTR/Daniel_results/SRS4245406.RData")
WT <- SRS4245406$WT
WT <- as.matrix(WT)
gKO = 'Cftr'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[which(g_wei == 0.1)]
length(gList)
