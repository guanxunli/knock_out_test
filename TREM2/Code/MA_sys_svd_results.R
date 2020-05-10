library(fgsea)
library(ggplot2)
library(enrichR)
source("NKX2-1/Code/utility.R")


#####################################################
## inpu data
X_1_both_0_T <- readRDS("NKX2-1/results/svd_results/MA_met1_both_0_TRUE.rds")
X_1_col_0_T <- readRDS("NKX2-1/results/svd_results/MA_met1_col_0_TRUE.rds")
X_1_row_0_T <- readRDS("NKX2-1/results/svd_results/MA_met1_row_0_TRUE.rds")
X_2_col_1_F <- readRDS("NKX2-1/results/svd_results/MA_met2_col_1_FALSE.rds")
X_2_col_1_T <- readRDS("NKX2-1/results/svd_results/MA_met2_col_1_TRUE.rds")
X_2_col_2_F <- readRDS("NKX2-1/results/svd_results/MA_met2_col_2_FALSE.rds")
X_2_col_2_T <- readRDS("NKX2-1/results/svd_results/MA_met2_col_2_TRUE.rds")
X_2_col_3_F <- readRDS("NKX2-1/results/svd_results/MA_met2_col_3_FALSE.rds")
X_2_col_3_T <- readRDS("NKX2-1/results/svd_results/MA_met2_col_3_TRUE.rds")
X_2_row_1_F <- readRDS("NKX2-1/results/svd_results/MA_met2_row_1_FALSE.rds")
X_2_row_1_T <- readRDS("NKX2-1/results/svd_results/MA_met2_row_1_TRUE.rds")
X_2_row_2_F <- readRDS("NKX2-1/results/svd_results/MA_met2_row_2_FALSE.rds")
X_2_row_2_T <- readRDS("NKX2-1/results/svd_results/MA_met2_row_2_TRUE.rds")
X_2_row_3_F <- readRDS("NKX2-1/results/svd_results/MA_met2_row_3_FALSE.rds")
X_2_row_3_T <- readRDS("NKX2-1/results/svd_results/MA_met2_row_3_TRUE.rds")

## marker genes
markerGenes <- read.csv('NKX2-1/Data/pnas.1906663116.sd01.csv', stringsAsFactors = FALSE, row.names = 1)
markerGenes$T.test.p.value <- as.numeric(markerGenes$T.test.p.value)
markerGenes <- markerGenes[complete.cases(markerGenes),]
markerGenesAT1 <- markerGenes$gene_short_name[markerGenes$T.test.p.value < 0.05]

markerGenes <- read.csv('NKX2-1/Data/pnas.1906663116.sd05.csv', stringsAsFactors = FALSE)
markerGenes <- markerGenes[,1:10]
markerGenes$T.test.p.value <- as.numeric(markerGenes$T.test.p.value)
markerGenes <- markerGenes[complete.cases(markerGenes),]
markerGenesAT2 <- markerGenes$gene_short_name[markerGenes$T.test.p.value < 0.05]

markerGenes <- unique(c(markerGenesAT1, markerGenesAT2)) # 7808
length(markerGenes)

######################################################################################################
################################## check results                 #####################################
######################################################################################################
out_2_row_1_F_all <- check_fun(X = X_2_row_1_F, d = 2, alpha = 0)
out_2_row_1_F_b <- check_fun(X = X_2_row_1_F, d = 2, alpha = 2)
out_2_col_1_F_all <- check_fun(X = X_2_col_1_F, d = 2, alpha = 0)
out_2_col_1_F_f <- check_fun(X = X_2_col_1_F, d = 2, alpha = 1)
