library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source("NKX2-1/Code/utility.R")

## inpu data
X_ima_row0 <- readRDS("NKX2-1/results/ima_results/Daniel_ima_row0.rds")
X_ima_col0 <- readRDS("NKX2-1/results/ima_results/Daniel_ima_col0.rds")
X_ima_both0 <- readRDS("NKX2-1/results/ima_results/Daniel_ima_both0.rds")

sum(Im(X_ima_row0)^2)
sum(Im(X_col_row0)^2)
sum(Im(X_both_row0)^2)

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

# dRegulation change

dRegulation <- function (manifoldOutput) {
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl("^X_", geneList)]
  geneList <- gsub("^X_", "", geneList)
  nGenes <- length(geneList)
  eGenes <- nrow(manifoldOutput)/2
  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl("^Y_", eGeneList)]
  eGeneList <- gsub("^Y_", "", eGeneList)
  if (nGenes != eGenes) {
    stop("Number of identified and expected genes are not the same")
  }
  if (!all(eGeneList == geneList)) {
    stop("Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order")
  }
  dMetric <- sapply(seq_len(nGenes), function(G) {
    X <- manifoldOutput[G, ]
    Y <- manifoldOutput[(G + nGenes), ]
    X_real <- Re(X)
    X_im <- Im(X)
    Y_real <- Re(Y)
    Y_im <- Im(Y)
    I_real <- rbind(X_real, Y_real)
    I_im  <- rbind(X_im, Y_im)
    O_real <- dist(I_real)
    O_im <- dist(I_im)
    O_real <- as.numeric(O_real)
    O_im <- as.numeric(O_im)
    O <- sqrt(O_real^2 + O_im^2)
    return(O)
  })
  lambdaValues <- seq(-2, 2, length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- MASS::boxcox(dMetric ~ 1, plot = FALSE, lambda = lambdaValues)
  BC <- BC$x[which.max(BC$y)]
  if (BC < 0) {
    nD <- 1/(dMetric^BC)
  }
  else {
    nD <- dMetric^BC
  }
  Z <- scale(nD)
  E <- mean(dMetric^2)
  FC <- dMetric^2/E
  pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = "fdr")
  dOut <- data.frame(gene = geneList, distance = dMetric, Z = Z, 
                     FC = FC, p.value = pValues, p.adj = pAdjusted)
  dOut <- dOut[order(dOut$p.value), ]
  dOut <- as.data.frame.array(dOut)
  return(dOut)
}


## fix P-values
# X is manifold alignment results
# d is the dimension of MA needed
# alpha = 0 means both 2d components, alpha = 1 means first, alpha = 2 means back d components
fixPValues <- function(X, d = 2, alpha = 1){
  if(alpha == 1){
    X <- X[, 1:d]
  } else if(alpha == 2){
    X <- X[, (ncol(X)/2 + 1):(ncol(X)/2 + d)]
  } else if(alpha == 0){
    X <- X[, c(1:d, (ncol(X)/2 + 1):(ncol(X)/2 + d))]
  }
  X <- list("manifoldAlignment" = X)
  X$diffRegulation <- dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-1]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
  return(X)
}

# check results for imaginary
out_row0 <- check_fun(X = X_ima_row0, d = 5, alpha = 1)


