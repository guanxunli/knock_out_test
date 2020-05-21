library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("HNF4A-HNF4G/Code/utility.R")

## inpu data
yan_met_1_0 <- readRDS("HNF4A-HNF4G/results/Yan_dir_results/yan_MA_row_Lt_1_0.rds")
yan_met_1_0 <- yan_met_1_0[, -1]
yan_met_1_1 <- readRDS("HNF4A-HNF4G/results/Yan_dir_results/yan_MA_row_Lt_1_1.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

# sensitivity
out_1_0 <- check_sensitivity(X = yan_met_1_0, d_index = c(1:10), alpha = 1)
saveRDS(out_1_0, "HNF4A-HNF4G/results/Yan_dir_results/out_1_0.rds")
out_1_1 <- check_sensitivity(X = yan_met_1_1, d_index = c(1:10), alpha = 1)
saveRDS(out_1_1, "HNF4A-HNF4G/results/Yan_dir_results/out_1_1.rds")


## plot sensitivity
out_1_0 <- readRDS("HNF4A-HNF4G/results/Yan_dir_results/out_1_0.rds")
plot_sencitivity(out_1_0)
out_1_1 <- readRDS("HNF4A-HNF4G/results/Yan_dir_results/out_1_1.rds")
plot_sencitivity(out_1_1)


## check with trivial method 

## row trivial
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
gList_row <- unique(c(gList_a, gList_g))

out_1_0 <- check_fun(X = yan_met_1_0, d = 2, alpha = 1)
check_intersect(gList_row, out_1_0$gene)

out_1_1 <- check_fun(X = yan_met_1_1, d = 2, alpha = 1)
check_intersect(gList_row, out_1_1$gene)

## column trivial
g_wei_a <- abs(WT[, 'Hnf4a'])
range(g_wei_a)
g_wei_g <- abs(WT[, 'Hnf4g'])
range(g_wei_g)
gList_a <- colnames(WT)[which(g_wei_a == 0.1)]
gList_g <- colnames(WT)[which(g_wei_g == 0.1)]
gList_col <- unique(c(gList_a, gList_g))

check_intersect(gList_col, out_1_0$gene)
check_intersect(gList_col, out_1_1$gene)

check_intersect(gList_col, out_1_0$gene)
check_intersect(gList_col, out_1_1$gene)

check_intersect(gList_row, gList_col)
check_intersect(out_1_0$gene, out_1_1$gene)
check_intersect(out_1_0$gene, out_ori_row0$gene)
check_intersect(out_1_1$gene, out_ori_row0$gene)

check_intersect(out_1_0$gene, out_row_sym$gene)



