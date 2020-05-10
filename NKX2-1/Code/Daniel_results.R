library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("NKX2-1/Code/utility.R")

## inpu data
X_ori_row0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_ori_row0.rds")
X_ori_col0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_ori_col0.rds")
X_ori_both0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_ori_both0.rds")
X_scc_row0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_scc_row0.rds")
X_scc_col0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_scc_col0.rds")
X_scc_both0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_scc_both0.rds")

X_ori_row0 <- X_ori_row0$MA
X_ori_col0 <- X_ori_col0$MA
X_ori_both0 <- X_ori_both0$MA
X_scc_row0 <- X_scc_row0$MA
X_scc_col0 <- X_scc_col0$MA
X_scc_both0 <- X_scc_both0$MA

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

out_ori_row0 <- check_fun(X = X_ori_row0, d = 2, alpha = 1)
out_ori_col0 <- check_fun(X = X_ori_col0, d = 2, alpha = 1)
