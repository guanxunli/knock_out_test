library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("TREM2/Code/utility.R")

## inpu data
yan_met_1_0 <- readRDS("TREM2/results/Yan_dir_results/Trem2yan_MA_row_Lt_1_0.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

# sensitivity
out_1_0 <- check_sensitivity(X = yan_met_1_0, d_index = c(1:10), alpha = 1)
saveRDS(out_1_0, "TREM2/results/Yan_dir_results/out_1_0.rds")


## plot sensitivity
out_1_0 <- readRDS("TREM2/results/Yan_dir_results/out_1_0.rds")
plot_sencitivity(out_1_0)

## check with trivial method 

## row trivial
load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'
g_wei <- abs(WT[gKO, ])
gList_row <- colnames(WT)[which(g_wei > 0.2)]

out_1_0 <- check_fun(X = yan_met_1_0, d = 2, alpha = 1)
check_intersect(gList_row, out_1_0$gene)


## column trivial
g_wei <- abs(WT[, gKO])
gList_col <- colnames(WT)[which(g_wei > 0.2)]
check_intersect(gList_col, out_1_0$gene)

check_intersect(gList_row, gList_col)
check_intersect(out_1_0$gene, out_ori_row0$gene)

check_intersect(out_1_0$gene, out_row_sym$gene)

