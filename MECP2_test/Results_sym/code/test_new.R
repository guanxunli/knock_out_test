library(Matrix)
library(scTenifoldNet)

manifoldAlignment_daniel1 <- function(X, Y, d = 10, lambda = 0.9){
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
  W <- -W
  diag(W) <- -apply(W, 2, sum)
  
  W <- (W + t(W)) / 2

  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-4]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

manifoldAlignment_daniel2 <- function(X, Y, d = 10, lambda = 0.9){
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
  W <- -W
  diag(W) <- -apply(W, 2, sum)
  
  W <- t(W) %*% W

  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-4]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

manifoldAlignment_daniel3 <- function(X, Y, d = 10, lambda = 0.9){
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
  W <- -W
  diag(W) <- -apply(W, 2, sum)
  
  W <- W %*% t(W)

  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-4]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

manifoldAlignment_yan1 <- function(X, Y, d = 5, lambda = 0.9){
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
  W <- Matrix(W)

  W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  W <- (W + t(W)) / 2

  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-4]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

manifoldAlignment_yan2 <- function(X, Y, d = 5, lambda = 0.9){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  n = dim(X)[1]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X + 1
  wY <- Y + 1
  wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  diag(W) <- 0

  Drow <- apply(abs(W), 1, sum)
  Dcol <- c(apply(abs(X),2,sum), apply(abs(Y),2,sum)) + mean(c(apply(abs(X),2,sum), apply(abs(Y),2,sum)))* lambda

  W <- Matrix(W)

  W = diag(Dcol) - W * (Drow^(-1)) * Dcol
  W = t(W) %*% W

  E <- suppressWarnings(RSpectra::eigs(W, d*3, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-4]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

# ## load the first data set
# load("SRS3059998.RData")

# WT <- SRS3059998$WT
# KO <- SRS3059998$KO

# set.seed(1)
# MA <- manifoldAlignment_daniel1(WT, KO, d = 5)
# saveRDS(MA, "daniel1_SRS3059998.rds")

# set.seed(1)
# MA <- manifoldAlignment_daniel2(WT, KO, d = 5)
# saveRDS(MA, "daniel2_SRS3059998.rds")

# set.seed(1)
# MA <- manifoldAlignment_daniel3(WT, KO, d = 5)
# saveRDS(MA, "daniel3_SRS3059998.rds")

# set.seed(1)
# MA <- manifoldAlignment_yan1(WT, KO, d = 5)
# saveRDS(MA, "yan1_SRS3059998.rds")

# set.seed(1)
# MA <- manifoldAlignment_yan2(WT, KO, d = 5)
# saveRDS(MA, "yan2_SRS3059998.rds")


## load the second data set
load("SRS3059999.RData")
WT <- SRS3059999$WT
KO <- SRS3059999$KO

set.seed(1)
MA <- manifoldAlignment_daniel1(WT, KO, d = 5)
saveRDS(MA, "daniel1_SRS3059999.rds")

set.seed(1)
MA <- manifoldAlignment_daniel2(WT, KO, d = 5)
saveRDS(MA, "daniel2_SRS3059999.rds")

set.seed(1)
MA <- manifoldAlignment_daniel3(WT, KO, d = 5)
saveRDS(MA, "daniel3_SRS3059999.rds")

set.seed(1)
MA <- manifoldAlignment_yan1(WT, KO, d = 5)
saveRDS(MA, "yan1_SRS3059999.rds")

set.seed(1)
MA <- manifoldAlignment_yan2(WT, KO, d = 5)
saveRDS(MA, "yan2_SRS3059999.rds")