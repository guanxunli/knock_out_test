library(ggplot2)
library(scTenifoldNet)
library(RSpectra)

## input data
load('NKX2-1/Daniel_Results/GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)

## define function
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
  E <- suppressWarnings(RSpectra::eigs(W, d * 2, "SM"))
  E$values <- suppressWarnings(Mod(E$values))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[, E$values > 1e-08]
  alignedNet <- E$vectors[, 1:d]
  colnames(alignedNet) <- paste0("NLMA ", seq_len(d))
  rownames(alignedNet) <- c(paste0("X_", sharedGenes), 
                            paste0("Y_", sharedGenes))
  return(alignedNet)
}

## row as 0
KO <- WT
gKO <- 'Nkx2-1'
KO[gKO,] <- 0
set.seed(1)
MA <- manifoldAlignment(WT, KO)
saveRDS(MA, "Daniel_ima_row0.rds")

## column as 0
KO <- WT
KO[, gKO] <- 0
set.seed(1)
MA <- manifoldAlignment(WT, KO)
saveRDS(MA, "Daniel_ima_col0.rds")

## both as 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
set.seed(1)
MA <- manifoldAlignment(WT, KO)
saveRDS(MA, "Daniel_ima_both0.rds")
