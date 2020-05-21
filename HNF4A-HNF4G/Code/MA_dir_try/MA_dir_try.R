library(fgsea)
library(ggplot2)
library(scTenifoldNet)
library(Matrix)

load("HNF4A-HNF4G/Daniel_Results/GSM3477499.RData")
WT <- GSM3477499$WT
WT <- as.matrix(WT)
gKO = c('Hnf4a','Hnf4g')


## modify manifold alignment
## normalized symmetry method
manifoldAlignment_norm_sym <- function(X, Y, d = 30, transpose = F){
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
  
  Wd = sqrt(rowSums(W))
  W =  t(t(W / Wd) / Wd)
  W <- as(W, "sparseMatrix")  
  
  E <- suppressWarnings(RSpectra::eigs_sym(W, d*2, which = "LA"))
  E$values <- suppressWarnings(as.numeric(E$values))
  neworder <- order(E$values, decreasing = TRUE)
  E$values <- E$values[neworder]
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[, neworder]
  E$vectors <- E$vectors[,E$values < (1-1e-8)]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

manifoldAlignment_nnorm_sym <- function(X, Y, d = 30, transpose = F){
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
  
  E <- suppressWarnings(RSpectra::eigs_sym(W, d*2, which = "SM"))
  E$values <- suppressWarnings(as.numeric(E$values))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
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
MA_row_nnorm_sym <- manifoldAlignment_nnorm_sym(X = WT, Y = KO, d = 50)
saveRDS(MA_row_nnorm_sym, "MA_row_nnorm_sym.rds")
# col 0 
KO <- WT
KO[, gKO] <- 0
MA_col_nnorm_sym <- manifoldAlignment_nnorm_sym(X = WT, Y = KO, d = 50)
saveRDS(MA_col_nnorm_sym, "MA_col_nnorm_sym.rds")
# both 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
MA_both_nnorm_sym <- manifoldAlignment_nnorm_sym(X = WT, Y = KO, d = 50)
saveRDS(MA_both_nnorm_sym, "MA_both_nnorm_sym.rds")

## normalized
# row 0 
KO <- WT
KO[gKO, ] <- 0
MA_row_norm_sym <- manifoldAlignment_norm_sym(X = WT, Y = KO, d = 50)
saveRDS(MA_row_norm_sym, "MA_row_norm_sym.rds")
# col 0 
KO <- WT
KO[, gKO] <- 0
MA_col_norm_sym <- manifoldAlignment_norm_sym(X = WT, Y = KO, d = 50)
saveRDS(MA_col_norm_sym, "MA_col_norm_sym.rds")
# both 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
MA_both_norm_sym <- manifoldAlignment_norm_sym(X = WT, Y = KO, d = 50)
saveRDS(MA_both_norm_sym, "MA_both_norm_sym.rds")

