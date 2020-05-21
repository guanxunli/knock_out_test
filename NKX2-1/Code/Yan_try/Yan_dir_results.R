library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("NKX2-1/Code/utility.R")

## inpu data
yan_met_1_0 <- readRDS("NKX2-1/results/Yan_dir_results/yan_MA_row_Lt_1_0.rds")
yan_met_1_1 <- readRDS("NKX2-1/results/Yan_dir_results/yan_MA_row_Lt_1_1.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

# sensitivity
out_1_0 <- check_sensitivity(X = yan_met_1_0, d_index = c(1:10), alpha = 1)
saveRDS(out_1_0, "NKX2-1/results/Yan_dir_results/out_1_0.rds")
out_1_1 <- check_sensitivity(X = yan_met_1_1, d_index = c(1:10), alpha = 1)
saveRDS(out_1_1, "NKX2-1/results/Yan_dir_results/out_1_1.rds")


## plot sensitivity
out_1_0 <- readRDS("NKX2-1/results/Yan_dir_results/out_1_0.rds")
plot_sencitivity(out_1_0)
out_1_1 <- readRDS("NKX2-1/results/Yan_dir_results/out_1_1.rds")
plot_sencitivity(out_1_1)


## check with trivial method 

## row trivial
load('NKX2-1/Daniel_Results/GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)
gKO <- 'Nkx2-1'
g_wei <- abs(WT[gKO, ])
gList_row <- colnames(WT)[which(g_wei > 0.4)]

out_1_0 <- check_fun(X = yan_met_1_0, d = 2, alpha = 1)
check_intersect(gList_row, out_1_0$gene)

out_1_1 <- check_fun(X = yan_met_1_1, d = 3, alpha = 1)
check_intersect(gList_row, out_1_1$gene)

## column trivial
g_wei <- abs(WT[, gKO ])
gList_col <- colnames(WT)[which(g_wei > 0.4)]
check_intersect(gList_col, out_1_0$gene)
check_intersect(gList_col, out_1_1$gene)

check_intersect(gList_row, gList_col)
check_intersect(out_1_0$gene, out_1_1$gene)
check_intersect(out_1_0$gene, out_ori_row0$gene)
check_intersect(out_1_1$gene, out_ori_row0$gene)

check_intersect(out_1_0$gene, out_row_sym$gene)



