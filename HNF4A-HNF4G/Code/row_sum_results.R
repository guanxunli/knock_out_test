library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("HNF4A-HNF4G/Code/utility.R")
## inpu data
MA_row_sym <- readRDS("HNF4A-HNF4G/results/row_sum_results//MA_row_sym.rds")
MA_col_sym <- readRDS("HNF4A-HNF4G/results/row_sum_results/MA_col_sym.rds")
MA_both_sym <- readRDS("HNF4A-HNF4G/results/row_sum_results/MA_both_sym.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

# sensitivity
out_row_sym <- check_sensitivity(X = MA_row_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_row_sym, "HNF4A-HNF4G/results/row_sum_results/out_row_sym.rds")
out_col_sym <- check_sensitivity(X = MA_col_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_col_sym, "HNF4A-HNF4G/results/row_sum_results/out_col_sym.rds")
out_both_sym <- check_sensitivity(X = MA_both_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_both_sym, "HNF4A-HNF4G/results/row_sum_results/out_both_sym.rds")

## plot sensitivity
out_row_sym <- readRDS("HNF4A-HNF4G/results/row_sum_results/out_row_sym.rds")
plot_sencitivity(out_row_sym)
out_col_sym <- readRDS("HNF4A-HNF4G/results/row_sum_results/out_col_sym.rds")
plot_sencitivity(out_col_sym)
out_both_sym <- readRDS("HNF4A-HNF4G/results/row_sum_results/out_both_sym.rds")
plot_sencitivity(out_both_sym)

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
length(gList)

out_row_sym <- check_fun(X = MA_row_sym, d = 2, alpha = 1)
n <- length(out_row_sym$gene)
print(c(n, length(intersect(gList, out_row_sym$gene))))

out_col_sym <- check_fun(X = MA_col_sym, d = 2, alpha = 1)
n <- length(out_col_sym$gene)
print(c(n, length(intersect(gList, out_col_sym$gene))))

out_both_sym <- check_fun(X = MA_both_sym, d = 2, alpha = 1)
n <- length(out_both_sym$gene)
print(c(n, length(intersect(gList, out_both_sym$gene))))

check_intersect(out_row_sym$gene, out_col_sym$gene)
check_intersect(out_row_sym$gene, out_both_sym$gene)
check_intersect(out_col_sym$gene, out_both_sym$gene)

check_intersect(out_row_sym$gene, out_ori_row0$gene)
check_intersect(out_row_sym$gene, out_ori_col0$gene)
check_intersect(out_row_sym$gene, out_ori_both0$gene)

check_intersect(out_col_sym$gene, out_ori_row0$gene)
check_intersect(out_col_sym$gene, out_ori_col0$gene)
check_intersect(out_col_sym$gene, out_ori_both0$gene)

check_intersect(out_both_sym$gene, out_ori_row0$gene)
check_intersect(out_both_sym$gene, out_ori_col0$gene)
check_intersect(out_both_sym$gene, out_ori_both0$gene)



