library(ggplot2)
library(scTenifoldNet)
library(RSpectra)

## input data
load('GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)
gKO <- 'Nkx2-1'

manifoldAlignment_yan <- function(X, Y, d = 30, weight_method = 1, lambda = 0.9, sqrt1 = 0){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  n = dim(X)[1]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X+1
  wY <- Y+1
  wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  diag(W) <- 0
  Drow <- apply(W, 1, sum)
  Dcol <- apply(W, 2, sum)
  W <- as(W, "sparseMatrix")  
  
  if(sqrt1 == 1){
    Dcol = Dcol^(0.5)
  }
  
  if(weight_method == 1){
    W = diag(Dcol) - W * (Drow^(-1)) * Dcol
    #  diag(Dcol) -diag(Dcol) %*% diag(Drow^(-1)) %*% W
  }else{
    Dcol[c(1:n)+n] = Dcol[c(1:n)]
    W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  }
  
  W = t(W) %*% W
  
  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
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

KO <- WT
KO[gKO, ] <- 0
weight_method <- 1
sqrt1 <- 0
MA_row_sym <- manifoldAlignment_yan(X = WT, Y = KO, d = 50, weight_method = weight_method, sqrt1 = sqrt1)
saveRDS(MA_row_sym, paste0(gKO,paste0("yan_MA_row_Lt_",weight_method,"_",sqrt1,".rds")))



