library(fgsea)
library(ggplot2)
library(scTenifoldNet)

## input data
load('GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)

# manifoldAlignment <- function(X, Y, d = 30){
#   sharedGenes <- intersect(rownames(X), rownames(Y))
#   X <- X[sharedGenes, sharedGenes]
#   Y <- Y[sharedGenes, sharedGenes]
#   X <- X+1
#   Y <- Y+1
#   # wX <-  X %*% t(X) #+t(X) %*% X +
#   # wY <- Y %*% t(Y) #+ t(Y) %*% Y +
#   # diag(wX) <- 0
#   # diag(wY) <- 0
#   # wX = wX / max(wX)
#   # wY = wY / max(wY)
#   alpha = 0.5
#   beta = 0.5
#   indegreeX = rowSums(X) + 1
#   outdegreeX = colSums(X) + 1
#   Do = diag(outdegreeX^(-alpha))
#   Di = diag(indegreeX^(-beta))
#   Bd = Do %*% X %*% Di %*% t(X) %*%Do
#   Cd = Di %*% t(X) %*% Do %*% X %*%Di
#   wX = Bd + Cd
#
#   indegreeY = rowSums(Y) + 1
#   outdegreeY = colSums(Y) + 1
#   Do = diag(outdegreeY^(-alpha))
#   Di = diag(indegreeY^(-beta))
#   Bd = Do %*% Y %*% Di %*% t(Y) %*%Do
#   Cd = Di %*% t(Y) %*% Do %*% Y %*%Di
#   wY = Bd + Cd
#
#   wX = wX / max(wX)
#   wY = wY / max(wY)
#   L <- diag(length(sharedGenes))
#
#   wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
#   W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
#   # W <- (W + t(W))/2
#   W <- -W
#   diag(W) <- 0
#   diag(W) <- -apply(W, 2, sum)
#   # E <- suppressWarnings(RSpectra::eigs_sym(W, d*2))
#   E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "SR"))
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

manifoldAlignment <- function(X, Y, d = 30){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  X <- X+1
  Y <- Y+1
  wX <- X %*% t(X) +t(X) %*% X
  wY <- Y %*% t(Y) + t(Y) %*% Y
  diag(wX) <- 0
  diag(wY) <- 0
  wX = wX / max(abs(wX))
  wY = wY / max(abs(wY))

  L <- diag(length(sharedGenes))

  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  # W <- (W + t(W))/2
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
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

#### setting row as 0
gKO <- 'Nkx2-1'
KO <- WT
KO[gKO, ] <- 0
MA_row0 <- manifoldAlignment(WT, KO, d = 50)
saveRDS(MA_row0, "Yan_row0.rds")
rm(MA_row0)

#### setting column as 0
gKO <- 'Nkx2-1'
KO <- WT
KO[, gKO] <- 0
MA_col0 <- manifoldAlignment(WT, KO, d = 50)
saveRDS(MA_col0, "Yan_col0.rds")
rm(MA_col0)

#### setting both row and column as 0
gKO <- 'Nkx2-1'
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
MA_both0 <- manifoldAlignment(WT, KO, d = 50)
saveRDS(MA_both0, "Yan_both0.rds")
rm(MA_both0)

#### both row and column divide by 2
gKO <- 'Nkx2-1'
KO <- WT
KO[, gKO] <- KO[, gKO]/2
KO[gKO, ] <- KO[gKO, ]/2
MA_both2 <- manifoldAlignment(WT, KO, d = 50)
saveRDS(MA_both2, "Yan_both2.rds")
rm(MA_both2)
