library(ggplot2)
library(scTenifoldNet)
library(RSpectra)

## input data
load('NKX2-1/Daniel_Results/GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)
gKO <- 'Nkx2-1'

## not normalized method
manifoldAlignment_nnorm <- function(X, Y, d = 30, transpose = F){
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

## normalized method
manifoldAlignment_norm <- function(X, Y, d = 30, transpose = F){
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

## not normalized
# row 0 
KO <- WT
KO[gKO, ] <- 0
MA_row_nnorm <- manifoldAlignment_nnorm(X = WT, Y = KO, d = 50)
saveRDS(MA_row_nnorm, "MA_row_nnorm.rds")
# col 0 
KO <- WT
KO[, gKO] <- 0
MA_row_nnorm <- manifoldAlignment_nnorm(X = WT, Y = KO, d = 50)
saveRDS(MA_row_nnorm, "MA_col_nnorm.rds")
# both 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
MA_row_nnorm <- manifoldAlignment_nnorm(X = WT, Y = KO, d = 50)
saveRDS(MA_row_nnorm, "MA_both_nnorm.rds")

## normalized
# row 0 
KO <- WT
KO[gKO, ] <- 0
MA_row_nnorm <- manifoldAlignment_norm(X = WT, Y = KO, d = 50)
saveRDS(MA_row_nnorm, "MA_row_norm.rds")
# col 0 
KO <- WT
KO[, gKO] <- 0
MA_row_nnorm <- manifoldAlignment_norm(X = WT, Y = KO, d = 50)
saveRDS(MA_row_nnorm, "MA_col_norm.rds")
# both 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
MA_row_nnorm <- manifoldAlignment_norm(X = WT, Y = KO, d = 50)
saveRDS(MA_row_nnorm, "MA_both_norm.rds")

