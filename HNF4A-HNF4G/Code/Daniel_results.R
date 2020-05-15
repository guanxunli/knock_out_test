library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("HNF4A-HNF4G/Code/utility.R")

## inpu data
X_ori_row0 <- readRDS("HNF4A-HNF4G/results/Daniel_results/Daniel_ori_row0.rds")
X_ori_col0 <- readRDS("HNF4A-HNF4G/results/Daniel_results/Daniel_ori_col0.rds")
X_ori_both0 <- readRDS("HNF4A-HNF4G/results/Daniel_results/Daniel_ori_both0.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

## check sensitivity
out_ori_row <- check_sensitivity(X = X_ori_row0, d_index = c(1:10))
saveRDS(out_ori_row, file = "HNF4A-HNF4G/results/Daniel_results/out_ori_row.rds")
out_ori_col <- check_sensitivity(X = X_ori_col0, d_index = c(1:10))
saveRDS(out_ori_col, file = "HNF4A-HNF4G/results/Daniel_results/out_ori_col.rds")
out_ori_both <- check_sensitivity(X = X_ori_both0, d_index = c(1:10))
saveRDS(out_ori_both, file = "HNF4A-HNF4G/results/Daniel_results/out_ori_both.rds")

## plot sensitivity
out_ori_row0 <- readRDS("HNF4A-HNF4G/results/Daniel_results/out_ori_row.rds")
plot_sencitivity(out_ori_row0)
out_ori_col0 <- readRDS("HNF4A-HNF4G/results/Daniel_results/out_ori_col.rds")
plot_sencitivity(out_ori_col0)
out_ori_both0 <- readRDS("HNF4A-HNF4G/results/Daniel_results/out_ori_both.rds")
plot_sencitivity(out_ori_both0)

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

out_ori_row0 <- check_fun(X = X_ori_row0, d = 2, alpha = 1)
n <- length(out_ori_row0$gene)
print(c(n, length(intersect(gList, out_ori_row0$gene))))

out_ori_col0 <- check_fun(X = X_ori_col0, d = 2, alpha = 1)
n <- length(out_ori_col0$gene)
print(c(n, length(intersect(gList, out_ori_col0$gene))))

out_ori_both0 <- check_fun(X = X_ori_both0, d = 2, alpha = 1)
n <- length(out_ori_both0$gene)
print(c(n, length(intersect(gList, out_ori_both0$gene))))

check_intersect(out_ori_row0$gene, out_ori_col0$gene)
check_intersect(out_ori_row0$gene, out_ori_both0$gene)
check_intersect(out_ori_col0$gene, out_ori_both0$gene)
