library(fgsea)
library(ggplot2)
library(Matrix)

load("HNF4A-HNF4G/Daniel_Results/GSM3477499.RData")
WT <- GSM3477499$WT
WT <- as.matrix(WT)
gKO = c('Hnf4a','Hnf4g')

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
  E <- suppressWarnings(RSpectra::eigs(W, d * 2, "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors, 2, as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[, E$values > 1e-08]
  alignedNet <- E$vectors[, seq_len(d)]
  colnames(alignedNet) <- paste0("NLMA ", seq_len(d))
  rownames(alignedNet) <- c(paste0("X_", sharedGenes), paste0("Y_", sharedGenes))
  return(alignedNet)
}

## row as 0
KO <- WT
KO[gKO,] <- 0
set.seed(1)
MA <- manifoldAlignment(WT, KO)
saveRDS(MA, "MA_row_sym.rds")

## column as 0
KO <- WT
KO[, gKO] <- 0
set.seed(1)
MA <- manifoldAlignment(WT, KO)
saveRDS(MA, "MA_col_sym.rds")

## both as 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
set.seed(1)
MA <- manifoldAlignment(WT, KO)
saveRDS(MA, "MA_both_sym.rds")
