library(ggplot2)
library(scTenifoldNet)
library(RSpectra)

## input data
load("SRS4245406.RData")
WT <- SRS4245406$WT
WT <- as.matrix(WT)
gKO = 'Cftr'

manifoldAlignment <- function(X, Y, d = 30){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  wX <- X + 1
  wY <- Y + 1
  L <- diag(length(sharedGenes))
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
  W <- (W + t(W))/2
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


## not normalized
# row 0 
KO <- WT
KO[gKO, ] <- 0
MA_row_sym <- manifoldAlignment(X = WT, Y = KO, d = 50)
saveRDS(MA_row_sym, "MA_row_sym.rds")

# col 0 
KO <- WT
KO[, gKO] <- 0
MA_col_sym <- manifoldAlignment(X = WT, Y = KO, d = 50)
saveRDS(MA_col_sym, "MA_col_sym.rds")
# both 0
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
MA_both_sym <- manifoldAlignment(X = WT, Y = KO, d = 50)
saveRDS(MA_both_sym, "MA_both_sym.rds")


