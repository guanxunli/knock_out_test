library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("NKX2-1/Code/utility.R")

## inpu data
MA_both_nnorm <- readRDS("NKX2-1/results/MA_dir_results/MA_both_nnorm.rds")
MA_row_nnorm <- readRDS("NKX2-1/results/MA_dir_results/MA_row_nnorm.rds")
MA_col_nnorm <- readRDS("NKX2-1/results/MA_dir_results/MA_col_nnorm.rds")
MA_both_norm <- readRDS("NKX2-1/results/MA_dir_results/MA_both_norm.rds")
MA_row_norm <- readRDS("NKX2-1/results/MA_dir_results/MA_row_norm.rds")
MA_col_norm <- readRDS("NKX2-1/results/MA_dir_results/MA_col_norm.rds")

# marker gene
markerGenes <- read.csv('NKX2-1/Data/pnas.1906663116.sd01.csv', stringsAsFactors = FALSE, row.names = 1)
markerGenes$T.test.p.value <- as.numeric(markerGenes$T.test.p.value)
markerGenes <- markerGenes[complete.cases(markerGenes),]
markerGenesAT1 <- markerGenes$gene_short_name[markerGenes$T.test.p.value < 0.05]

markerGenes <- read.csv('NKX2-1//Data/pnas.1906663116.sd05.csv', stringsAsFactors = FALSE)
markerGenes <- markerGenes[,1:10]
markerGenes$T.test.p.value <- as.numeric(markerGenes$T.test.p.value)
markerGenes <- markerGenes[complete.cases(markerGenes),]
markerGenesAT2 <- markerGenes$gene_short_name[markerGenes$T.test.p.value < 0.05]

markerGenes <- unique(c(markerGenesAT1, markerGenesAT2)) # 7808
length(markerGenes)

######################################################################################################
################################## check results                 #####################################
######################################################################################################

out_row_nnorm <- check_fun(X = MA_row_nnorm, d = 5, alpha = 1)
out_col_nnorm <- check_fun(X = MA_col_nnorm, d = 5, alpha = 2)
out_both_nnorm <- check_fun(X = MA_both_nnorm, d = 5, alpha = 0)

out_row_norm <- check_fun(X = MA_row_norm, d = 2, alpha = 2)
out_col_norm <- check_fun(X = MA_col_norm, d = 5, alpha = 1)
out_both_norm <- check_fun(X = MA_both_norm, d = 5, alpha = 0)
