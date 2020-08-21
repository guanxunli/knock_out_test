library(scTenifoldNet)
library(Matrix)

SERGIO <- read.csv('simulation data/SERGIO_create_sim_data/simulationOutput.csv', header = FALSE)
rownames(SERGIO) <- paste0('G', seq_len(nrow(SERGIO)))
colnames(SERGIO) <- paste0('C', seq_len(ncol(SERGIO)))

countMatrix <- SERGIO
set.seed(1)
X <- scTenifoldNet::makeNetworks(countMatrix, q = 0.9)
set.seed(1)
X <- scTenifoldNet::tensorDecomposition(X)
X <- X$X
X <- as.matrix(X)

## check function
check_intersect <- function(gList1, gList2){
  return(c(length(gList1), length(gList2), length(intersect(gList1, gList2))))
}

# ## daniel
# manifoldAlignment <- function(X, Y, d = 10, lambda = 0.9){
#   sharedGenes <- intersect(rownames(X), rownames(Y))
#   X <- X[sharedGenes, sharedGenes]
#   n = dim(X)[1]
#   Y <- Y[sharedGenes, sharedGenes]
#   L <- diag(length(sharedGenes))
#   wX <- X+1
#   wY <- Y+1
#   wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
#   W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
#   diag(W) <- 0
#   W <- -W
#   diag(W) <- -apply(W, 2, sum)
#   
#   # W <- (W + t(W)) / 2
#   # W <- W %*% t(W)
#   # W <- t(W) %*% W
# 
#   E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
#   E$values <- suppressWarnings(as.numeric(E$values))
#   E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
#   newOrder <- order(E$values)
#   E$values <- E$values[newOrder]
#   E$vectors <- E$vectors[,newOrder]
#   E$vectors <- E$vectors[,E$values > 1e-8]
#   alignedNet <- E$vectors[,seq_len(d)]
#   colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
#   rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
#   return(alignedNet)
# }

# ## Yan method
# manifoldAlignment <- function(X, Y, d = 5, lambda = 0.9){
#   sharedGenes <- intersect(rownames(X), rownames(Y))
#   X <- X[sharedGenes, sharedGenes]
#   n = dim(X)[1]
#   Y <- Y[sharedGenes, sharedGenes]
#   L <- diag(length(sharedGenes))
#   wX <- X+1
#   wY <- Y+1
#   wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
#   W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
#   diag(W) <- 0
# 
#   Drow <- apply(W, 1, sum)
#   Dcol <- apply(W, 2, sum)
#   W <- Matrix(W)
# 
#   W = diag(Dcol) - W * (Drow^(-1)) * Dcol
#   # W = t(W) %*% W
#   W <- (W + t(W)) / 2
# 
#   E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
#   E$values <- suppressWarnings(as.numeric(E$values))
#   E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
#   newOrder <- order(E$values)
#   E$values <- E$values[newOrder]
#   E$vectors <- E$vectors[,newOrder]
#   E$vectors <- E$vectors[,E$values > 1e-8]
#   alignedNet <- E$vectors[,seq_len(d)]
#   colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
#   rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
#   return(alignedNet)
# }

# ## modify yan
# manifoldAlignment <- function(X, Y, d = 5, lambda = 0.9){
#   sharedGenes <- intersect(rownames(X), rownames(Y))
#   X <- X[sharedGenes, sharedGenes]
#   n = dim(X)[1]
#   Y <- Y[sharedGenes, sharedGenes]
#   L <- diag(length(sharedGenes))
#   wX <- X + 1
#   wY <- Y + 1
#   wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
#   W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
#   diag(W) <- 0
# 
#   Drow <- apply(abs(W), 1, sum)
#   # Dcol <- apply(abs(W), 2, sum)
#   Dcol <- c(apply(abs(X),2,sum), apply(abs(Y),2,sum)) + mean(c(apply(abs(X),2,sum), apply(abs(Y),2,sum)))* lambda
# 
#   W <- Matrix(W)
# 
#   W = diag(Dcol) - W * (Drow^(-1)) * Dcol
#   W = t(W) %*% W
# 
#   E <- suppressWarnings(RSpectra::eigs(W, d*3, 'SR'))
#   E$values <- suppressWarnings(as.numeric(E$values))
#   E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
#   newOrder <- order(E$values)
#   E$values <- E$values[newOrder]
#   E$vectors <- E$vectors[,newOrder]
#   E$vectors <- E$vectors[,E$values > 1e-8]
#   alignedNet <- E$vectors[,seq_len(d)]
#   colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
#   rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
#   return(alignedNet)
# }

## modify yan
dirtect_L <- function(x){
  x <- t(t(x) / colSums(x))
  e <- suppressWarnings(RSpectra::eigs(t(x), 1, 'LR'))
  Phi <- as.numeric(e$vectors)
  if (Phi[1] < 0){
    Phi <- - Phi
  }
  tmp1 <- diag(sqrt(Phi)) %*% x %*% diag(1/(sqrt(Phi)))
  tmp2 <- diag(1/sqrt(Phi)) %*% t(x) %*% diag(sqrt(Phi))
  L <- diag(1, nrow(x)) - (tmp1 + tmp2)/2
  return(L)
}

manifoldAlignment <- function(X, Y, d = 5, lambda = 1.5){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  n = dim(X)[1]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X + 1
  wY <- Y + 1
  
  
  ## different laplace
  # Lx <- dirtect_L(wX)
  # Ly <- dirtect_L(wY)
  # 
  # Lx_diag <- diag(Lx)
  # Ly_diag <- diag(Ly)
  # diag(Lx) <- 0
  # diag(Ly) <- 0
  # 
  # LXY <- lambda * (sum(Lx) + sum(Ly)) / (2 * sum(L)) * L
  # 
  # L <- rbind(cbind(Lx, LXY), cbind(t(LXY), Ly))
  # diag(L) <- -apply(L, 2, sum)
  
  ## same laplace
  wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  L <- dirtect_L(W)

  E <- suppressWarnings(RSpectra::eigs(L, d*2, 'SR'))
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

p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- manifoldAlignment(X,Y, d = 5)
  MA <- MA[, 1:3]
  DR <- dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

## test method
g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
check_intersect(g_list, g_true)

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
check_intersect(g_list, g_true)

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
check_intersect(g_list, g_true)

