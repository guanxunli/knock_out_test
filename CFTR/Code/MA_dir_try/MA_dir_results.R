library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("CFTR/Code/utility.R")

## inpu data
MA_both_nnorm <- readRDS("CFTR/results/MA_dir_results/MA_both_nnorm.rds")
MA_row_nnorm <- readRDS("CFTR/results/MA_dir_results/MA_row_nnorm.rds")
MA_col_nnorm <- readRDS("CFTR/results/MA_dir_results/MA_col_nnorm.rds")
MA_both_norm <- readRDS("CFTR/results/MA_dir_results/MA_both_norm.rds")
MA_row_norm <- readRDS("CFTR/results/MA_dir_results/MA_row_norm.rds")
MA_col_norm <- readRDS("CFTR/results/MA_dir_results/MA_col_norm.rds")

######################################################################################################
################################## check results                 #####################################
######################################################################################################


## sensitivity
out_row_nnorm <- check_sensitivity(X = MA_row_nnorm, d_index = c(1:10), alpha = 2)
saveRDS(out_row_nnorm, "CFTR/results/MA_dir_results/out_row_nnorm.rds")
out_col_nnorm <- check_sensitivity(X = MA_col_nnorm, d_index = c(1:10), alpha = 1)
saveRDS(out_col_nnorm, "CFTR/results/MA_dir_results/out_col_nnorm.rds")
out_both_nnorm <- check_sensitivity(X = MA_both_nnorm, d_index = c(1:10), alpha = 0)
saveRDS(out_both_nnorm, "CFTR/results/MA_dir_results/out_both_nnorm.rds")

## using p-value rather than p-adj
out_row_norm <- check_sensitivity(X = MA_row_norm, d_index = c(1:10), alpha = 2)
saveRDS(out_row_norm, "CFTR/results/MA_dir_results/out_row_norm.rds")
out_col_norm <- check_sensitivity(X = MA_col_norm, d_index = c(1:10), alpha = 1)
saveRDS(out_col_norm, "CFTR/results/MA_dir_results/out_col_norm.rds")
out_both_norm <- check_sensitivity(X = MA_both_norm, d_index = c(1:10), alpha = 0)
saveRDS(out_both_norm, "CFTR/results/MA_dir_results/out_both_norm.rds")

## plot results
out_row_norm <- readRDS("CFTR/results/MA_dir_results/out_row_norm.rds")
plot_sencitivity(out_row_norm)
out_col_norm <- readRDS("CFTR/results/MA_dir_results/out_col_norm.rds")
plot_sencitivity(out_col_norm)
out_both_norm <- readRDS("CFTR/results/MA_dir_results/out_both_norm.rds")
plot_sencitivity(out_both_norm)
out_row_nnorm <- readRDS("CFTR/results/MA_dir_results/out_row_nnorm.rds")
plot_sencitivity(out_row_nnorm)
out_col_nnorm <- readRDS("CFTR/results/MA_dir_results/out_col_nnorm.rds")
plot_sencitivity(out_col_nnorm)
out_both_nnorm <- readRDS("CFTR/results/MA_dir_results/out_both_nnorm.rds")
plot_sencitivity(out_both_nnorm)

## check first d and last d
out_row_nnorm_f <- check_fun(X = MA_row_nnorm, d = 2, alpha = 1)
length(out_row_nnorm_f$gene)
head(out_row_nnorm_f$gene)
out_row_nnorm_b <- check_fun(X = MA_row_nnorm, d = 2, alpha = 2)
head(out_row_nnorm_b$gene)
length(out_row_nnorm_b$gene)

out_col_nnorm_f <- check_fun(X = MA_col_nnorm, d = 2, alpha = 1)
length(out_col_nnorm_f$gene)
head(out_col_nnorm_f$gene)
out_col_nnorm_b <- check_fun(X = MA_col_nnorm, d = 2, alpha = 2)
head(out_col_nnorm_b$gene)
length(out_col_nnorm_b$gene)

out_both_nnorm_f <- check_fun(X = MA_both_nnorm, d = 2, alpha = 1)
length(out_both_nnorm_f$gene)
head(out_both_nnorm_f$gene)
out_both_nnorm_b <- check_fun(X = MA_both_nnorm, d = 2, alpha = 2)
head(out_both_nnorm_b$gene)
length(out_both_nnorm_b$gene)

out_row_norm_f <- check_fun(X = MA_row_norm, d = 2, alpha = 1)
length(out_row_norm_f$gene)
head(out_row_norm_f$gene)
out_row_norm_b <- check_fun(X = MA_row_norm, d = 2, alpha = 2)
head(out_row_norm_b$gene)
length(out_row_norm_b$gene)

out_col_norm_f <- check_fun(X = MA_col_norm, d = 2, alpha = 1)
length(out_col_norm_f$gene)
head(out_col_norm_f$gene)
out_col_norm_b <- check_fun(X = MA_col_norm, d = 2, alpha = 2)
head(out_col_norm_b$gene)
length(out_col_norm_b$gene)

out_both_norm_f <- check_fun(X = MA_both_norm, d = 2, alpha = 1)
length(out_both_norm_f$gene)
head(out_both_norm_f$gene)
out_both_norm_b <- check_fun(X = MA_col_norm, d = 2, alpha = 2)
head(out_both_norm_b$gene)
length(out_both_norm_b$gene)

## check with trivial method
load("CFTR/Daniel_results/SRS4245406.RData")
WT <- SRS4245406$WT
WT <- as.matrix(WT)
gKO = 'Cftr'
g_wei <- abs(WT[gKO, ])
gList <- colnames(WT)[order(g_wei, decreasing = TRUE)]

out_row_norm <- check_fun(X = MA_row_norm, d = 2, alpha = 2)
n <- length(out_row_norm$gene)
print(c(n, length(intersect(gList[1:n], out_row_norm$gene))))

out_col_norm <- check_fun(X = MA_col_norm, d = 2, alpha = 1)
n <- length(out_col_norm$gene)
print(c(n, length(intersect(gList[1:n], out_col_norm$gene))))

out_both_norm <- check_fun(X = MA_both_norm, d = 2, alpha = 2)
n <- length(out_both_norm$gene)
print(c(n, length(intersect(gList[1:n], out_both_norm$gene))))

out_row_nnorm <- check_fun(X = MA_row_nnorm, d = 5, alpha = 2)
n <- length(out_row_nnorm$gene)
print(c(n, length(intersect(gList[1:n], out_row_nnorm$gene))))

out_col_nnorm <- check_fun(X = MA_col_nnorm, d = 5, alpha = 1)
n <- length(out_col_nnorm$gene)
print(c(n, length(intersect(gList[1:n], out_col_nnorm$gene))))

out_both_nnorm <- check_fun(X = MA_both_nnorm, d = 5, alpha = 0)
n <- length(out_both_nnorm$gene)
print(c(n, length(intersect(gList[1:n], out_both_nnorm$gene))))

## check intersection
check_intersect(out_row_nnorm$gene, out_col_nnorm$gene)
check_intersect(out_row_nnorm$gene, out_both_nnorm$gene)
check_intersect(out_col_nnorm$gene, out_both_nnorm$gene)
check_intersect(out_row_nnorm$gene, out_ori_row0$gene)
check_intersect(out_row_nnorm$gene, out_ori_col0$gene)
check_intersect(out_row_nnorm$gene, out_ori_both0$gene)
check_intersect(out_col_nnorm$gene, out_ori_both0$gene)
check_intersect(out_col_nnorm$gene, out_ori_col0$gene)
check_intersect(out_col_nnorm$gene, out_ori_row0$gene)
check_intersect(out_both_nnorm$gene, out_ori_both0$gene)
check_intersect(out_both_nnorm$gene, out_ori_col0$gene)
check_intersect(out_both_nnorm$gene, out_ori_row0$gene)

