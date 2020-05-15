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

manifoldAlignment <- function(X, Y, d = 30){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  wX <- X + 1
  wY <- Y + 1
  L <- diag(length(sharedGenes))
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  # W <- (W + t(W))/2
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
  W <- (W + t(W))/2
  # E <- suppressWarnings(RSpectra::eigs_sym(W, d*2))
  E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

#### original method
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- scTenifoldNet::manifoldAlignment(X,Y)
  MA <- MA[, 1:3]
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
  MA <- manifoldAlignment(X,Y, d = 3)
  DR <- scTenifoldNet::dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
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

#### setting column as 0
p_value_KO <- function(gKO){
  X_sym <- (X + t(X))/2
  Y <- X_sym
  # Y <- X
  Y[, gKO] <- 0
  Y <- (Y + t(Y))/2
  MA <- manifoldAlignment(X,Y, d = 2)
  # MA <- scTenifoldNet::manifoldAlignment(X_sym,Y_sym, d = 10)
  DR <- scTenifoldNet::dRegulation(MA)
  # DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$FC <- DR$distance^2/mean(DR$distance^2)
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

#### setting both 0
p_value_KO <- function(gKO){
  X_sym <- (X + t(X))/2
  Y <- X_sym
  # Y <- X
  Y[, gKO] <- 0
  Y[gKO, ] <- 0
  Y <- (Y + t(Y))/2
  MA <- manifoldAlignment(X,Y, d = 10)
  # MA <- scTenifoldNet::manifoldAlignment(X_sym,Y_sym, d = 10)
  DR <- scTenifoldNet::dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
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

#### divide both row and column by 2
p_value_KO <- function(gKO){
  X_sym <- (X + t(X))/2
  Y <- X_sym
  # Y <- X
  Y[, gKO] <- Y[, gKO]/2
  Y[gKO, ] <- Y[gKO, ]/2
  Y <- (Y + t(Y))/2
  MA <- manifoldAlignment(X,Y, d = 10)
  # MA <- scTenifoldNet::manifoldAlignment(X_sym,Y_sym, d = 10)
  DR <- scTenifoldNet::dRegulation(MA, minFC = 0)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
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

g_res <- p_value_KO(41)
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

#### divide both row and column by 2 and check order
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- Y[, gKO]/2
  Y[gKO, ] <- Y[gKO, ]/2
  MA <- manifoldAlignment(X,Y, d = 10)
  DR <- scTenifoldNet::dRegulation(MA, minFC = 0)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  # DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene[1:25]
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(41)
g_list <- g_res$gene[1:40]
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene[1:20]
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))
