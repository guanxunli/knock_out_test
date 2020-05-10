library(fgsea)
library(ggplot2)
library(scTenifoldNet)
library(Matrix)

## input data
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

WT <- readMM('GSM3716703_Nkx2-1_control_scRNAseq_matrix.mtx.gz')
rownames(WT) <- read.csv('GSM3716703_Nkx2-1_control_scRNAseq_genes.tsv.gz', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]
colnames(WT) <- readLines('GSM3716703_Nkx2-1_control_scRNAseq_barcodes.tsv.gz')
WT <- scQC(WT)
WT <- WT[rowMeans(WT != 0) > 0.05,]
WT <- WT[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT), ignore.case = TRUE),]

sccNet <- function(X, q = 0.95, nCell = 500, nNet = 25, K = 2, denoiseNet = TRUE){
  nGenes <- nrow(X)
  gList <- rownames(X)
  oNet <- pbapply::pbsapply(seq_len(nNet), function(Z){
    tNet <- Matrix::t(X[,sample(seq_len(ncol(X)), nCell, replace = TRUE)])
    tNet <- cor(as.matrix(tNet), method = 'sp')
    while(any(is.na(tNet))){
      tNet <- Matrix::t(X[,sample(seq_len(ncol(X)), nCell, replace = TRUE)])
      tNet <- cor(as.matrix(tNet), method = 'sp')
    }
    tNet <- round(tNet,1)
    diag(tNet) <- 0
    tNet[abs(tNet) < quantile(abs(tNet), q, na.rm = TRUE)] <- 0
    tNet <- Matrix::Matrix(tNet)
    return(tNet)
  })

  aNet <- matrix(0, nGenes, nGenes)
  rownames(aNet) <- colnames(aNet) <- gList
  for(i in seq_along(oNet)){
    tNet <- oNet[[i]]
    if(denoiseNet){
      if(K == ncol(X)){
        tNet <- svd(tNet)
        tNet <- tNet$u %*% diag(tNet$d) %*% t(tNet$v)
        rownames(tNet) <- colnames(tNet) <- gList
      } else if(K == 1){
        tNet <- RSpectra::svds(tNet,K)
        tNet <- tNet$u %*% (tNet$d) %*% t(tNet$v)
        rownames(tNet) <- colnames(tNet) <- gList
      } else {
        tNet <- RSpectra::svds(tNet,K, maxitr = 1e6)
        tNet <- tNet$u %*% diag(tNet$d) %*% t(tNet$v)
        rownames(tNet) <- colnames(tNet) <- gList
      }

    }
    tNet <- Matrix::Matrix(tNet)
    aNet <- aNet + tNet[gList, gList]
  }
  aNet <- aNet/length(oNet)
  aNet <- Matrix::Matrix(aNet)
  return(aNet)
}

WT <- sccNet(WT, K = 10)
saveRDS(WT, file = 'net_try.rds')

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
