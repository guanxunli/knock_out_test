library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("CFTR/Code/utility.R")

## inpu data
MA_ori_row0 <- readRDS("CFTR/results/Daniel_results/Daniel_ori_row0.rds")
MA_ori_col0 <- readRDS("CFTR/results/Daniel_results/Daniel_ori_col0.rds")
MA_ori_both0 <- readRDS("CFTR/results/Daniel_results/Daniel_ori_both0.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

## check sensitivity
out_ori_row <- check_sensitivity(X = X_ori_row0, d_index = c(1:10))
saveRDS(out_ori_row, file = "TREM2/results/Daniel_results/out_ori_row.rds")
out_ori_col <- check_sensitivity(X = X_ori_col0, d_index = c(1:10))
saveRDS(out_ori_col, file = "TREM2/results/Daniel_results/out_ori_col.rds")
out_ori_both <- check_sensitivity(X = X_ori_both0, d_index = c(1:10))
saveRDS(out_ori_both, file = "TREM2/results/Daniel_results/out_ori_both.rds")

## plot sensitivity
out_ori_row0 <- readRDS("TREM2/results/Daniel_results/out_ori_row.rds")
plot_sencitivity(out_ori_row0)
out_ori_col0 <- readRDS("TREM2/results/Daniel_results/out_ori_col.rds")
plot_sencitivity(out_ori_col0)
out_ori_both0 <- readRDS("TREM2/results/Daniel_results/out_ori_both.rds")
plot_sencitivity(out_ori_both0)

## check head gene and trivial method
load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[order(g_wei, decreasing = TRUE)]

out_ori_row0 <- check_fun(X = X_ori_row0, d = 2, alpha = 1)
length(out_ori_row0$gene)
head(out_ori_row0$gene)
n <- length(out_ori_row0$gene)
print(c(n, length(intersect(gList[1:n], out_ori_row0$gene))))

out_ori_col0 <- check_fun(X = X_ori_col0, d = 2, alpha = 1)
length(out_ori_col0$gene)
head(out_ori_col0$gene)
n <- length(out_ori_col0$gene)
print(c(n, length(intersect(gList[1:n], out_ori_col0$gene))))

out_ori_both0 <- check_fun(X = X_ori_both0, d = 2, alpha = 1)
length(out_ori_both0$gene)
head(out_ori_both0$gene)
n <- length(out_ori_both0$gene)
print(c(n, length(intersect(gList[1:n], out_ori_both0$gene))))