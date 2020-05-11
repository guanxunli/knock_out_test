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

## method 1
manifoldAlignment1 <- function(X, Y, d = 30){
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  L <- diag(length(sharedGenes))
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  Wr <- rowSums(W)
  Wc <- colSums(W)
  Wnow <- diag(Wr^(-0.5)) %*% W %*% diag(Wc^(-0.5))
  E <- svds(Wnow, d * 2)
  select.index = E$d < 1- 10^(-8)
  alignedNet <- cbind(as.matrix(E$u[,select.index])[,seq_len(d)], as.matrix(E$v[,select.index])[,seq_len(d)])
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

## method 2
manifoldAlignment2 <- function(X, Y, d = 30, transpose = F){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)))
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)))
  
  L <- diag(length(sharedGenes )* 2)
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  Wd = rowSums(W)
  W <- diag(Wd) - W
  W <- as(W, "sparseMatrix")  
  
  E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

## method 3
manifoldAlignment_normalize <- function(X, Y, d = 30, transpose = F){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)))
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)))
  
  L <- diag(length(sharedGenes )* 2)
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  Wd = rowSums(W)
  W = diag(Wd^(-0.5)) %*% W %*% diag(Wd^(-0.5))
  W <- as(W, "sparseMatrix")  
  
  E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "LR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[,E$values < (1-1e-8)]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}


#############################################################try 1###############################################################################
### setting row as zero
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- manifoldAlignment1(X,Y, d = 2)
  MA <- MA[, 3:4]
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

#### setting column as 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  MA <- manifoldAlignment1(X,Y, d = 2)
  MA <- MA[, 1:2]
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

#### setting both 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  Y[gKO, ] <- 0
  MA <- manifoldAlignment1(X,Y, d = 2)
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

##################################################################try 2 #####################################################################
### setting row as zero
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  d <- 2
  MA <- manifoldAlignment2(X,Y, d = d, transpose = F)
  MA <- MA[, (d + 1): (2 * d)]
  # MA <- MA[, 1:d]
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

#### setting column as 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  MA <- manifoldAlignment2(X,Y, d = 2, normalize = 1)
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

#### setting both 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  Y[gKO, ] <- 0
  MA <- manifoldAlignment2(X,Y, d = 2, normalize = 0)
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

##################################################################try 3 #####################################################################
### setting row as zero
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  d <- 2
  MA <- manifoldAlignment_normalize(X,Y, d = d, transpose = F)
  MA <- MA[, (d + 1): (2 * d)]
  # MA <- MA[, 1:d]
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

#### setting column as 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  MA <- manifoldAlignment2(X,Y, d = 2, normalize = 1)
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

#### setting both 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  Y[gKO, ] <- 0
  MA <- manifoldAlignment2(X,Y, d = 2, normalize = 0)
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
