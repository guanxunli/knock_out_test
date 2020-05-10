library(ggplot2)
library(scTenifoldNet)
library(RSpectra)

## input data
load('NKX2-1/Daniel_Results/GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)

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

## normalize can be 0, 1,2  1 and 2 are two different way to normalize. 0 is not to normalize.
manifoldAlignment2 <- function(X, Y, d = 30, normalize = 0, transpose = F){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)) )
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)) )
  
  L <- diag(length(sharedGenes )* 2)
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  if(normalize == 1){
    Wd = rowSums(W)
    W = diag(Wd^(-0.5)) %*% W %*% diag(Wd^(-0.5))
    W <- -W
    diag(W) <- 1
  }else if(normalize == 2){
    Wd = rowSums(W)
    W = diag(Wd^(-1)) %*% W
    W <- -W
    diag(W) <- 1
  }else if(normalize == 3){
    Wd = colSums(W)
    W = W %*% diag(Wd^(-1))
    W <- -W
    diag(W) <- 1
  }else{
    W <- -W
    diag(W) <- 0
    diag(W) <- -apply(W, 2, sum)
  }
  
  # E <- suppressWarnings(RSpectra::eigs_sym(W, d*2))
  E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}


#### server try
# nohup Rscript  test.R 3 5 >  re.txt &
Args <- commandArgs()
met <- as.numeric(Args[6])
row_col <- as.numeric(Args[7])
norm_met <- as.numeric(Args[8])
trans_not <- as.numeric(Args[9])

gKO <- 'Nkx2-1'
KO <- WT

if (row_col == "both"){
  KO[gKO, ] <- 0
  KO[, gKO] <- 0
} else if(row_col == "row"){
  KO[gKO, ] <- 0
} else if (row_col == "col"){
  KO[, gKO] <- 0
}

if (met == 1){
  MA <- manifoldAlignment1(WT, KO, d = 50)
} else if (met == 2){
  MA <- manifoldAlignment2(WT, KO, d = 50, normalize = norm_met, transpose = trans_not)
}

saveRDS(MA, file = paste0("MA_met", met, "_", row_col, "_", norm_met, "_", trans_not, ".rds"))



