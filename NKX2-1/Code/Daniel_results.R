library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("NKX2-1/Code/utility.R")

## inpu data
X_ori_row0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_ori_row0.rds")
X_ori_col0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_ori_col0.rds")
X_ori_both0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_ori_both0.rds")
# X_scc_row0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_scc_row0.rds")
# X_scc_col0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_scc_col0.rds")
# X_scc_both0 <- readRDS("NKX2-1/results/Daniel_results/Daniel_scc_both0.rds")

X_ori_row0 <- X_ori_row0$MA
X_ori_col0 <- X_ori_col0$MA
X_ori_both0 <- X_ori_both0$MA
# X_scc_row0 <- X_scc_row0$MA
# X_scc_col0 <- X_scc_col0$MA
# X_scc_both0 <- X_scc_both0$MA

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

# sensitivity
out_ori_row0 <- check_sensitivity(X = X_ori_row0, d_index = c(1:10), alpha = 1)
saveRDS(out_ori_row0, "NKX2-1/results/Daniel_results/out_ori_row.rds")
out_ori_col0 <- check_sensitivity(X = X_ori_col0, d_index = c(1:10), alpha = 1)
saveRDS(out_ori_col0, "NKX2-1/results/Daniel_results/out_ori_col.rds")
out_ori_both0 <- check_sensitivity(X = X_ori_both0, d_index = c(1:10), alpha = 1)
saveRDS(out_ori_both0, "NKX2-1/results/Daniel_results/out_ori_both.rds")

## plot sensitivity
out_ori_row0 <- readRDS("NKX2-1/results/Daniel_results/out_ori_row.rds")
plot_sencitivity(out_ori_row0)
out_ori_col0 <- readRDS("NKX2-1/results/Daniel_results/out_ori_col.rds")
plot_sencitivity(out_ori_col0)
out_ori_both0 <- readRDS("NKX2-1/results/Daniel_results/out_ori_both.rds")
plot_sencitivity(out_ori_both0)

## check head gene and trivial method
load('NKX2-1/Daniel_Results/GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)
gKO <- 'Nkx2-1'
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

check_intersect(out_ori_row0$gene, out_ori_col0$gene)
check_intersect(out_ori_row0$gene, out_ori_both0$gene)
check_intersect(out_ori_col0$gene, out_ori_both0$gene)


