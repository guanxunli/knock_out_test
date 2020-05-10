library(fgsea)
library(ggplot2)
library(scTenifoldNet)

## input data
load('GSM3716703.RData')
WT <- GSM3716703$WT
WT <- as.matrix(WT)

manifoldAlignment <- function(X, Y, d = 30){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X+1
  wY <- Y+1
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  W <- (W + t(W))/2
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
  E <- suppressWarnings(RSpectra::eigs_sym(W, d*2))
  # E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "SR"))
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
MA_row0 <- manifoldAlignment(WT, KO)
DR_row0 <- dRegulation(MA_row0, minFC = 0)
saveRDS(MA_row0, "MA_row0.rds")
saveRDS(DR_row0, "DR_row0.rds")
rm(MA_row0, DR_row0)

#### setting column as 0
gKO <- 'Nkx2-1'
KO <- WT
KO[, gKO] <- 0
MA_col0 <- manifoldAlignment(WT, KO)
DR_col0 <- dRegulation(MA_col0, minFC = 0)
saveRDS(MA_col0, "MA_col0.rds")
saveRDS(DR_col0, "DR_col0.rds")
rm(MA_col0, DR_col0)

#### setting both row and column as 0
gKO <- 'Nkx2-1'
KO <- WT
KO[gKO, ] <- 0
KO[, gKO] <- 0
MA_both0 <- manifoldAlignment(WT, KO)
DR_both0 <- dRegulation(MA_both0, minFC = 0)
saveRDS(MA_both0, "MA_both0.rds")
saveRDS(DR_both0, "DR_both0.rds")
rm(MA_both0, DR_both0)

#### both row and column divide by 2
gKO <- 'Nkx2-1'
KO <- WT
KO[, gKO] <- KO[, gKO]/2
KO[gKO, ] <- KO[gKO, ]/2
MA_both2 <- manifoldAlignment(WT, KO)
DR_both2 <- dRegulation(MA_both2, minFC = 0)
saveRDS(MA_both2, "MA_both2.rds")
saveRDS(DR_both2, "DR_both2.rds")
rm(MA_both2, DR_both2)
