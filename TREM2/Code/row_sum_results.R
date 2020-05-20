library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("TREM2/Code/utility.R")
## inpu data
MA_row_sym <- readRDS("TREM2/results/row_sum_results//MA_row_sym.rds")
MA_col_sym <- readRDS("TREM2/results/row_sum_results/MA_col_sym.rds")
MA_both_sym <- readRDS("TREM2/results/row_sum_results/MA_both_sym.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

# sensitivity
out_ori_row0 <- check_sensitivity(X = MA_row_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_ori_row0, "TREM2/results/row_sum_results/out_ori_row.rds")
out_ori_col0 <- check_sensitivity(X = MA_col_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_ori_col0, "TREM2/results/row_sum_results/out_ori_col.rds")
out_ori_both0 <- check_sensitivity(X = MA_both_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_ori_both0, "TREM2/results/row_sum_results/out_ori_both.rds")

## plot sensitivity
out_row_sym <- readRDS("TREM2/results/row_sum_results/out_ori_row.rds")
plot_sencitivity(out_row_sym)
out_col_sym <- readRDS("TREM2/results/row_sum_results/out_ori_col.rds")
plot_sencitivity(out_col_sym)
out_both_sym <- readRDS("TREM2/results/row_sum_results/out_ori_both.rds")
plot_sencitivity(out_both_sym)

## check with trivial method
load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[which(g_wei > 0.2)]

## check head gene and trivial method
load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'
g_wei <- abs(WT[, gKO ])
gList <- colnames(WT)[which(g_wei > 0.2)]

out_row_sym2 <- check_fun(X = MA_row_sym, d = 2, alpha = 1)
length(out_row_sym$gene)
head(out_row_sym$gene)
n <- length(out_row_sym$gene)
print(c(n, length(intersect(gList, out_row_sym$gene))))

out_col_sym2 <- check_fun(X = MA_col_sym, d = 2, alpha = 1)
length(out_col_sym$gene)
head(out_col_sym$gene)
n <- length(out_col_sym$gene)
print(c(n, length(intersect(gList, out_col_sym$gene))))

out_both_sym2 <- check_fun(X = MA_both_sym, d = 2, alpha = 1)
length(out_both_sym$gene)
head(out_both_sym$gene)
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



