library(ggplot2)
library(statsExpressions)
library(patchwork)
library(circlize)
library(pbapply)

SERGIO <- read.csv('SERGIO/simulationOutput.csv', header = FALSE)
rownames(SERGIO) <- paste0('G', seq_len(nrow(SERGIO)))
colnames(SERGIO) <- paste0('C', seq_len(ncol(SERGIO)))

countMatrix <- SERGIO
set.seed(1)

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

X <- sccNet(countMatrix, K = 5)
X <- as.matrix(X)

manifoldAlignment <- function(X, Y, d = 30){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- cov(X)
  wY <- cov(Y)
  wX <- wX/max(abs(wX))
  wY <- wY/max(abs(wY))
  wX <- wX + 1
  wY <- wY + 1
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  W <- (W + t(W))/2
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 1, sum)
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

#### original method
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- scTenifoldNet::manifoldAlignment(X,Y, d = 10)
  DR <- scTenifoldNet::dRegulation(MA, minFC = 0)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

#### setting row as zero
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- manifoldAlignment(X,Y, d = 10)
  DR <- scTenifoldNet::dRegulation(MA, minFC = 0)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

#### setting column as 0
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- 0
  MA <- manifoldAlignment(X,Y, d = 10)
  DR <- scTenifoldNet::dRegulation(MA, minFC = 0)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

#### setting both 0
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  Y[, gKO] <- 0
  MA <- manifoldAlignment(X,Y, d = 10)
  DR <- scTenifoldNet::dRegulation(MA, minFC = 0)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

#### divide both row and column by 2
p_value_KO <- function(gKO){
  Y <- X
  Y[, gKO] <- Y[, gKO]/2
  Y[gKO, ] <- Y[gKO, ]/2
  MA <- manifoldAlignment(X,Y, d = 10)
  DR <- scTenifoldNet::dRegulation(MA, minFC = 0)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  # DR$FC <- DR$distance^2/mean(DR$distance^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(41)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
length(g_list)
length(g_true)
length(intersect(g_list, g_true))
