## check with trivial method
load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[order(g_wei, decreasing = TRUE)]
source("TREM2/Code/utility.R")

# trivial method
g_number <- seq(10, 150, by = 10)
n <- length(g_number)
enrich_res <- rep(0, n)
for (i in 1:n){
  gList_use <- gList[1:g_number[i]]
  E <- enrichFunction(gList_use, fdrThreshold = 0.05, nCategories = 20)
  enrich_res[i] <- length(unique(E$Term))
}
saveRDS(enrich_res, "TREM2/results/trivial_res.rds")
library(ggplot2)
enrich_res <- readRDS("TREM2/results/trivial_res.rds")
ggplot(data = NULL, aes(x = g_number, y = enrich_res)) + geom_line() + geom_point()
