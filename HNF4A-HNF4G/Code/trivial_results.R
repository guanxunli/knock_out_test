## check head gene and trivial method
load("HNF4A-HNF4G/Daniel_Results/GSM3477499.RData")
WT <- GSM3477499$WT
WT <- as.matrix(WT)
gKO = c('Hnf4a','Hnf4g')
g_wei_a <- abs(WT['Hnf4a', ])
range(g_wei_a)
g_wei_g <- abs(WT['Hnf4g', ])
range(g_wei_g)
gList_a <- colnames(WT)[which(g_wei_a == 0.1)]
gList_g <- colnames(WT)[which(g_wei_g == 0.1)]
gList <- unique(c(gList_a, gList_g))
source("HNF4A-HNF4G/Code/utility.R")
E <- enrichFunction(gList, fdrThreshold = 0.05, nCategories = 20)
length(unique(E$Term))

