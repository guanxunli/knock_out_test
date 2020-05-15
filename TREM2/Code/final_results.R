library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("TREM2/Code/utility.R")

## input data
MA_both_nnorm_sym <- readRDS("TREM2/results/MA_dir_results/MA_both_nnorm_sym.rds")
MA_row_nnorm_sym <- readRDS("TREM2/results/MA_dir_results/MA_row_nnorm_sym.rds")
MA_col_nnorm_sym <- readRDS("TREM2/results/MA_dir_results/MA_col_nnorm_sym.rds")
MA_both_norm_sym <- readRDS("TREM2/results/MA_dir_results/MA_both_norm_sym.rds")
MA_row_norm_sym <- readRDS("TREM2/results/MA_dir_results/MA_row_norm_sym.rds")
MA_col_norm_sym <- readRDS("TREM2/results/MA_dir_results/MA_col_norm_sym.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################

## sensitivity
out_row_nnorm_sym <- check_sensitivity_fb(X = MA_row_nnorm_sym, d_index = c(1:10)) 
saveRDS(out_row_nnorm_sym, "TREM2/results/Final_results/out_row_nnorm_sym.rds")
out_col_nnorm_sym <- check_sensitivity_fb(X = MA_col_nnorm_sym, d_index = c(1:10))
saveRDS(out_col_nnorm_sym, "TREM2/results/Final_results/out_col_nnorm_sym.rds")
out_both_nnorm_sym <- check_sensitivity_fb(X = MA_both_nnorm_sym, d_index = c(1:10))
saveRDS(out_both_nnorm_sym, "TREM2/results/Final_results/out_both_nnorm_sym.rds")

out_row_norm_sym <- check_sensitivity_fb(X = MA_row_norm_sym, d_index = c(1:10))
saveRDS(out_row_norm_sym, "TREM2/results/Final_results/out_row_norm_sym.rds")
out_col_norm_sym <- check_sensitivity_fb(X = MA_col_norm_sym, d_index = c(1:10))
saveRDS(out_col_norm_sym, "TREM2/results/Final_results/out_col_norm_sym.rds")
out_both_norm_sym <- check_sensitivity_fb(X = MA_both_norm_sym, d_index = c(1:10))
saveRDS(out_both_norm_sym, "TREM2/results/Final_results/out_both_norm_sym.rds")

## plot results
out_row_norm_sym <- readRDS("TREM2/results/Final_results/out_row_norm_sym.rds")
plot_sencitivity(out_row_norm_sym)
out_col_norm_sym <- readRDS("TREM2/results/Final_results/out_col_norm_sym.rds")
plot_sencitivity(out_col_norm_sym)
out_both_norm_sym <- readRDS("TREM2/results/Final_results/out_both_norm_sym.rds")
plot_sencitivity(out_both_norm_sym)
out_row_nnorm_sym <- readRDS("TREM2/results/Final_results/out_row_nnorm_sym.rds")
plot_sencitivity(out_row_nnorm_sym)
out_col_nnorm_sym <- readRDS("TREM2/results/Final_results/out_col_nnorm_sym.rds")
plot_sencitivity(out_col_nnorm_sym)
out_both_nnorm_sym <- readRDS("TREM2/results/Final_results/out_both_nnorm_sym.rds")
plot_sencitivity(out_both_nnorm_sym)

## check with trivial method
load("TREM2/Daniel_results/GSE130626.RData")
WT <- GSE130626$WT
WT <- as.matrix(WT)
gKO <- 'Trem2'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[which(g_wei > 0.2)]

out_row_norm_sym <- check_fun_fb(X = MA_row_norm_sym, d = 2)
n <- length(out_row_norm_sym$gene)
print(c(n, length(intersect(gList, out_row_norm_sym$gene))))

out_col_norm_sym <- check_fun_fb(X = MA_col_norm_sym, d = 2)
n <- length(out_col_norm_sym$gene)
print(c(n, length(intersect(gList, out_col_norm_sym$gene))))

out_both_norm_sym <- check_fun_fb(X = MA_both_norm_sym, d = 2)
n <- length(out_both_norm_sym$gene)
print(c(n, length(intersect(gList, out_both_norm_sym$gene))))

out_row_nnorm_sym <- check_fun_fb(X = MA_row_nnorm_sym, d = 2)
n <- length(out_row_nnorm_sym$gene)
print(c(n, length(intersect(gList, out_row_nnorm_sym$gene))))

out_col_nnorm_sym <- check_fun_fb(X = MA_col_nnorm_sym, d = 2)
n <- length(out_col_nnorm_sym$gene)
print(c(n, length(intersect(gList, out_col_nnorm_sym$gene))))

out_both_nnorm_sym <- check_fun_fb(X = MA_both_nnorm_sym, d = 2)
n <- length(out_both_nnorm_sym$gene)
print(c(n, length(intersect(gList, out_both_nnorm_sym$gene))))

## check intersection
check_intersect(out_row_norm_sym$gene, out_col_norm_sym$gene)
check_intersect(out_row_norm_sym$gene, out_both_norm_sym$gene)
check_intersect(out_col_norm_sym$gene, out_both_norm_sym$gene)

check_intersect(out_row_nnorm_sym$gene, out_col_nnorm_sym$gene)
check_intersect(out_row_nnorm_sym$gene, out_both_nnorm_sym$gene)
check_intersect(out_col_nnorm_sym$gene, out_both_nnorm_sym$gene)

check_intersect(out_row_nnorm_sym$gene, out_row_norm_sym$gene)
check_intersect(out_row_nnorm_sym$gene, out_col_norm_sym$gene)
check_intersect(out_row_nnorm_sym$gene, out_both_norm_sym$gene)

check_intersect(out_col_nnorm_sym$gene, out_row_norm_sym$gene)
check_intersect(out_col_nnorm_sym$gene, out_col_norm_sym$gene)
check_intersect(out_col_nnorm_sym$gene, out_both_norm_sym$gene)

check_intersect(out_both_nnorm_sym$gene, out_row_norm_sym$gene)
check_intersect(out_both_nnorm_sym$gene, out_col_norm_sym$gene)
check_intersect(out_both_nnorm_sym$gene, out_both_norm_sym$gene)

check_intersect(out_row_nnorm_sym$gene, out_ori_row0$gene)
check_intersect(out_row_nnorm_sym$gene, out_ori_col0$gene)
check_intersect(out_row_nnorm_sym$gene, out_ori_both0$gene)

check_intersect(out_col_nnorm_sym$gene, out_ori_row0$gene)
check_intersect(out_col_nnorm_sym$gene, out_ori_col0$gene)
check_intersect(out_col_nnorm_sym$gene, out_ori_both0$gene)

check_intersect(out_both_nnorm_sym$gene, out_ori_row0$gene)
check_intersect(out_both_nnorm_sym$gene, out_ori_col0$gene)
check_intersect(out_both_nnorm_sym$gene, out_ori_both0$gene)

check_intersect(out_both_norm_sym$gene, out_ori_row0$gene)
check_intersect(out_both_norm_sym$gene, out_ori_col0$gene)
check_intersect(out_both_norm_sym$gene, out_ori_both0$gene)

check_intersect(out_col_norm_sym$gene, out_ori_row0$gene)
check_intersect(out_col_norm_sym$gene, out_ori_col0$gene)
check_intersect(out_col_norm_sym$gene, out_ori_both0$gene)

check_intersect(out_both_norm_sym$gene, out_ori_row0$gene)
check_intersect(out_both_norm_sym$gene, out_ori_col0$gene)
check_intersect(out_both_norm_sym$gene, out_ori_both0$gene)


