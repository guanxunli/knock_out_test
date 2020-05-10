library(ggplot2)
library(scTenifoldNet)

SERGIO <- read.csv('simulation data/SERGIO_create_sim_data/simulationOutput.csv', header = FALSE)
rownames(SERGIO) <- paste0('G', seq_len(nrow(SERGIO)))
colnames(SERGIO) <- paste0('C', seq_len(ncol(SERGIO)))

countMatrix <- SERGIO
set.seed(1)
X <- scTenifoldNet::makeNetworks(countMatrix, q = 0.9)
X <- scTenifoldNet::tensorDecomposition(X)
X <- X$X
X <- as.matrix(X)

manifoldAlignment <- function (X, Y, d = 30) {
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X + 1
  wY <- Y + 1
  wXY <- 0.9 * (sum(wX) + sum(wY))/(2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
  E <- suppressWarnings(RSpectra::eigs(W, d * 2, "SM"))
  E$values <- suppressWarnings(Mod(E$values))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[, E$values > 1e-08]
  alignedNet <- E$vectors[, 1:d]
  colnames(alignedNet) <- paste0("NLMA ", seq_len(d))
  rownames(alignedNet) <- c(paste0("X_", sharedGenes), 
                            paste0("Y_", sharedGenes))
  return(alignedNet)
}

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

#### original method
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- scTenifoldNet::manifoldAlignment(X,Y, d = 30)
  MA <- MA[, 1:10]
  DR <- scTenifoldNet::dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

#### setting row as zero
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- manifoldAlignment(X,Y, d = 10)
  DR <- dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

#### original setting column as 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  MA <- scTenifoldNet::manifoldAlignment(X, Y)
  MA <- MA[, 1:10]
  DR <- scTenifoldNet::dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

#### setting column 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  MA <- manifoldAlignment(X, Y)
  MA <- MA[, 1:10]
  DR <- dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))
