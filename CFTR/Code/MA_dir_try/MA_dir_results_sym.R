library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("CFTR/Code/utility.R")

## inpu data
MA_both_nnorm_sym <- readRDS("CFTR/results/MA_dir_results/MA_both_nnorm_sym.rds")
MA_row_nnorm_sym <- readRDS("CFTR/results/MA_dir_results/MA_row_nnorm_sym.rds")
MA_col_nnorm_sym <- readRDS("CFTR/results/MA_dir_results/MA_col_nnorm_sym.rds")
MA_both_norm_sym <- readRDS("CFTR/results/MA_dir_results/MA_both_norm_sym.rds")
MA_row_norm_sym <- readRDS("CFTR/results/MA_dir_results/MA_row_norm_sym.rds")
MA_col_norm_sym <- readRDS("CFTR/results/MA_dir_results/MA_col_norm_sym.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################


## sensitivity
out_row_nnorm_sym <- check_sensitivity(X = MA_row_nnorm_sym, d_index = c(1:10), alpha = 2)
saveRDS(out_row_nnorm_sym, "CFTR/results/MA_dir_results/out_row_nnorm_sym.rds")
out_col_nnorm_sym <- check_sensitivity(X = MA_col_nnorm_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_col_nnorm_sym, "CFTR/results/MA_dir_results/out_col_nnorm_sym.rds")
out_both_nnorm_sym <- check_sensitivity(X = MA_both_nnorm_sym, d_index = c(1:10), alpha = 0)
saveRDS(out_both_nnorm_sym, "CFTR/results/MA_dir_results/out_both_nnorm_sym.rds")

out_row_norm_sym <- check_sensitivity(X = MA_row_norm_sym, d_index = c(1:10), alpha = 2)
saveRDS(out_row_norm_sym, "CFTR/results/MA_dir_results/out_row_norm_sym.rds")
out_col_norm_sym <- check_sensitivity(X = MA_col_norm_sym, d_index = c(1:10), alpha = 1)
saveRDS(out_col_norm_sym, "CFTR/results/MA_dir_results/out_col_norm_sym.rds")
out_both_norm_sym <- check_sensitivity(X = MA_both_norm_sym, d_index = c(1:10), alpha = 0)
saveRDS(out_both_norm_sym, "CFTR/results/MA_dir_results/out_both_norm_sym.rds")

## plot results
out_row_norm_sym <- readRDS("CFTR/results/MA_dir_results/out_row_norm_sym.rds")
plot_sencitivity(out_row_norm_sym)
out_col_norm_sym <- readRDS("CFTR/results/MA_dir_results/out_col_norm_sym.rds")
plot_sencitivity(out_col_norm_sym)
out_both_norm_sym <- readRDS("CFTR/results/MA_dir_results/out_both_norm_sym.rds")
plot_sencitivity(out_both_norm_sym)
out_row_nnorm_sym <- readRDS("CFTR/results/MA_dir_results/out_row_nnorm_sym.rds")
plot_sencitivity(out_row_nnorm_sym)
out_col_nnorm_sym <- readRDS("CFTR/results/MA_dir_results/out_col_nnorm_sym.rds")
plot_sencitivity(out_col_nnorm_sym)
out_both_nnorm_sym <- readRDS("CFTR/results/MA_dir_results/out_both_nnorm_sym.rds")
plot_sencitivity(out_both_nnorm_sym)

## check first d and last d
out_row_nnorm_f_sym <- check_fun(X = MA_row_nnorm_sym, d = 2, alpha = 1)
length(out_row_nnorm_f_sym$gene)
head(out_row_nnorm_f_sym$gene)
out_row_nnorm_b_sym <- check_fun(X = MA_row_nnorm_sym, d = 2, alpha = 2)
head(out_row_nnorm_b_sym$gene)
length(out_row_nnorm_b_sym$gene)

out_col_nnorm_f_sym <- check_fun(X = MA_col_nnorm_sym, d = 2, alpha = 1)
length(out_col_nnorm_f_sym$gene)
head(out_col_nnorm_f_sym$gene)
out_col_nnorm_b_sym <- check_fun(X = MA_col_nnorm_sym, d = 2, alpha = 2)
head(out_col_nnorm_b_sym$gene)
length(out_col_nnorm_b_sym$gene)

out_both_nnorm_f_sym <- check_fun(X = MA_both_nnorm_sym, d = 2, alpha = 1)
length(out_both_nnorm_f_sym$gene)
head(out_both_nnorm_f_sym$gene)
out_both_nnorm_b_sym <- check_fun(X = MA_both_nnorm_sym, d = 2, alpha = 2)
head(out_both_nnorm_b_sym$gene)
length(out_both_nnorm_b_sym$gene)

out_row_norm_f_sym <- check_fun(X = MA_row_norm_sym, d = 2, alpha = 1)
length(out_row_norm_f_sym$gene)
head(out_row_norm_f_sym$gene)
out_row_norm_b_sym <- check_fun(X = MA_row_norm_sym, d = 2, alpha = 2)
head(out_row_norm_b_sym$gene)
length(out_row_norm_b_sym$gene)

out_col_norm_f_sym <- check_fun(X = MA_col_norm_sym, d = 2, alpha = 1)
length(out_col_norm_f_sym$gene)
head(out_col_norm_f_sym$gene)
out_col_norm_b_sym <- check_fun(X = MA_col_norm_sym, d = 2, alpha = 2)
head(out_col_norm_b_sym$gene)
length(out_col_norm_b_sym$gene)

out_both_norm_f_sym <- check_fun(X = MA_both_norm_sym, d = 2, alpha = 1)
length(out_both_norm_f_sym$gene)
head(out_both_norm_f_sym$gene)
out_both_norm_b_sym <- check_fun(X = MA_col_norm_sym, d = 2, alpha = 2)
head(out_both_norm_b_sym$gene)
length(out_both_norm_b_sym$gene)

## check with trivial method
load("CFTR/Daniel_results/SRS4245406.RData")
WT <- SRS4245406$WT
WT <- as.matrix(WT)
gKO = 'Cftr'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[which(g_wei == 0.1)]
length(gList)


out_row_norm_sym <- check_fun(X = MA_row_norm_sym, d = 2, alpha = 2)
n <- length(out_row_norm_sym$gene)
print(c(n, length(intersect(gList, out_row_norm_sym$gene))))

out_col_norm_sym <- check_fun(X = MA_col_norm_sym, d = 2, alpha = 1)
n <- length(out_col_norm_sym$gene)
print(c(n, length(intersect(gList, out_col_norm_sym$gene))))

out_both_norm_sym <- check_fun(X = MA_both_norm_sym, d = 2, alpha = 0)
n <- length(out_both_norm_sym$gene)
print(c(n, length(intersect(gList, out_both_norm_sym$gene))))

out_row_nnorm_sym <- check_fun(X = MA_row_nnorm_sym, d = 2, alpha = 2)
n <- length(out_row_nnorm_sym$gene)
print(c(n, length(intersect(gList, out_row_nnorm_sym$gene))))

out_col_nnorm_sym<- check_fun(X = MA_col_nnorm_sym, d = 2, alpha = 1)
n <- length(out_col_nnorm_sym$gene)
print(c(n, length(intersect(gList, out_col_nnorm_sym$gene))))

out_both_nnorm_sym <- check_fun(X = MA_both_nnorm_sym, d = 2, alpha = 0)
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


