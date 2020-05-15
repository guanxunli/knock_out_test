## modify dRegulation
dRegulation <- function (manifoldOutput) 
{
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl("^X_", geneList)]
  geneList <- gsub("^X_", "", geneList)
  nGenes <- length(geneList)
  eGenes <- nrow(manifoldOutput)/2
  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl("^Y_", eGeneList)]
  eGeneList <- gsub("^Y_", "", eGeneList)
  if (nGenes != eGenes) {
    stop("Number of identified and expected genes are not the same")
  }
  if (!all(eGeneList == geneList)) {
    stop("Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order")
  }
  dMetric <- sapply(seq_len(nGenes), function(G) {
    X <- manifoldOutput[G, ]
    Y <- manifoldOutput[(G + nGenes), ]
    I <- rbind(X, Y)
    O <- dist(I)
    O <- as.numeric(O)
    if (O == 0){
      O = 1e-17
    }
    return(O)
  })
  lambdaValues <- seq(-2, 2, length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- MASS::boxcox(dMetric ~ 1, plot = FALSE, lambda = lambdaValues)
  BC <- BC$x[which.max(BC$y)]
  if (BC < 0) {
    nD <- 1/(dMetric^BC)
  }
  else {
    nD <- dMetric^BC
  }
  Z <- scale(nD)
  E <- mean(dMetric^2)
  FC <- dMetric^2/E
  pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = "fdr")
  dOut <- data.frame(gene = geneList, distance = dMetric, Z = Z, 
                     FC = FC, p.value = pValues, p.adj = pAdjusted)
  dOut <- dOut[order(dOut$p.value), ]
  dOut <- as.data.frame.array(dOut)
  return(dOut)
}

## modify manifold alignment
## normalized method
manifoldAlignment_norm <- function(X, Y, d = 30, transpose = F){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)))
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)))
  
  L <- diag(length(sharedGenes )* 2)
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  # wXY <- 0.5 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  Wd = sqrt(rowSums(W))
  W =  t(t(W / Wd) / Wd)
  W <- as(W, "sparseMatrix")  
  
  E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "LR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  neworder <- order(E$values, decreasing = TRUE)
  E$values <- E$values[neworder]
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[, neworder]
  E$vectors <- E$vectors[,E$values < (1-1e-8)]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

## not normalized method
manifoldAlignment_nnorm <- function(X, Y, d = 30, transpose = F){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)))
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)))
  
  L <- diag(length(sharedGenes )* 2)
  # wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  wXY <- 0.5 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  Wd = rowSums(W)
  W <- diag(Wd) - W
  W <- as(W, "sparseMatrix")  
  
  E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

## normalized symmetry method
manifoldAlignment_norm_sym <- function(X, Y, d = 30, transpose = F){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)))
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)))
  
  L <- diag(length(sharedGenes )* 2)
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  Wd = sqrt(rowSums(W))
  W =  t(t(W / Wd) / Wd)
  W <- as(W, "sparseMatrix")  
  
  E <- suppressWarnings(RSpectra::eigs_sym(W, d*2, which = "LA"))
  E$values <- suppressWarnings(as.numeric(E$values))
  neworder <- order(E$values, decreasing = TRUE)
  E$values <- E$values[neworder]
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[, neworder]
  E$vectors <- E$vectors[,E$values < (1-1e-8)]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

manifoldAlignment_nnorm_sym <- function(X, Y, d = 30, transpose = F){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)))
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)))
  
  L <- diag(length(sharedGenes )* 2)
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  Wd = rowSums(W)
  W <- diag(Wd) - W
  W <- as(W, "sparseMatrix")  
  
  E <- suppressWarnings(RSpectra::eigs_sym(W, d*2, which = "SM"))
  E$values <- suppressWarnings(as.numeric(E$values))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

manifoldAlignment_new <- function(X, Y, d = 30, transpose = F, lambda = 0.9, normalize1 = 0){
  library(Matrix)
  library(RSpectra)
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  
  if(transpose){
    X = t(X)
    Y = t(Y)
  }
  
  wX <- X+1
  wY <- Y+1
  diag(wX) = 0
  diag(wY) = 0
  
  if(normalize1 == 1){
    rX = rowSums(wX)^(-1/2)
    cX = colSums(wX)^(-1/2)
    wX = diag(rX) %*% wX %*% diag(cX)
    rY = rowSums(wY)^(-1/2)
    cY = colSums(wY)^(-1/2)
    wY = diag(rY) %*% wY %*% diag(cY)
  }
  
  nsharedGenes = length(sharedGenes)
  wX = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wX), cbind(t(wX), matrix(0, nsharedGenes,nsharedGenes)))
  wY = rbind(cbind(matrix(0, nsharedGenes,nsharedGenes), wY), cbind(t(wY), matrix(0, nsharedGenes,nsharedGenes)))
  
  L <- diag(length(sharedGenes )* 2)
  wXY <- lambda * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  
  Wd = rowSums(W)
  W <- diag(Wd) - W
  W <- as(W, "sparseMatrix")  
  
  # E <- suppressWarnings(RSpectra::eigs(W, d*2, which = "SR"))
  E <- suppressWarnings(RSpectra::eigs_sym(W, d*2, which = "SM"))
  
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  
  alignedNet <- rbind( cbind(E$vectors[1:nsharedGenes,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes,seq_len(d)]),
                       cbind(E$vectors[(1:nsharedGenes) + nsharedGenes * 2,seq_len(d)], E$vectors[(1:nsharedGenes) + nsharedGenes * 3,seq_len(d)]))
  colnames(alignedNet) <- c(paste0('NLMA_u_ ', seq_len(d)),paste0('NLMA_v_ ', seq_len(d)))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}


## check results
# X is network
# gKO is knock out gene
# alpha to decide which part

fix_pvalue <- function(X, gKO, d, alpha = 1, beta = 1, normalize = FALSE){
  Y <- X
  
  if (beta == 1){
    Y[gKO, ] <- 0
  } else if(beta == 2){
    Y[, gKO] <- 0
  } else{
    Y[gKO, ]<- 0
    Y[, gKO] <- 0
  }

  if (normalize == TRUE){
    MA <- manifoldAlignment_norm(X, Y, d = d)
  } else{
    MA <- manifoldAlignment_nnorm(X,Y, d = d)
  }

  if (alpha == 1){
    MA <- MA[, 1:d]
  } else if (alpha == 2) {
    MA <- MA[, (d + 1): (2 * d)]
  }
  DR <- dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

fix_pvalue_sym <- function(X, gKO, d, alpha = 1, beta = 1, normalize = FALSE){
  Y <- X
  
  if (beta == 1){
    Y[gKO, ] <- 0
  } else if(beta == 2){
    Y[, gKO] <- 0
  } else{
    Y[gKO, ]<- 0
    Y[, gKO] <- 0
  }
  
  if (alpha == 1){
    MA <- MA[, 1:d]
  } else if(alpha == 2){
    MA <- MA[, (d + 1): (2 * d)]
  } 
  
  DR <- dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

fix_pvalue_new <- function(X, gKO, d = 30, alpha = 1, beta = 1, normalize = 1, lambda = 0.9){
  Y <- X
  
  if (beta == 1){
    Y[gKO, ] <- 0
  } else if(beta == 2){
    Y[, gKO] <- 0
  } else{
    Y[gKO, ]<- 0
    Y[, gKO] <- 0
  }
  
  MA <- manifoldAlignment_new(X, Y, d = d, transpose = F, lambda = lambda, normalize1 = normalize)

  
  if (alpha == 1){
    MA <- MA[, 1:d]
  } else if (alpha == 2) {
    MA <- MA[, (d + 1): (2 * d)]
  }
  DR <- dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}


## check function
check_intersect <- function(gList1, gList2){
  return(c(length(gList1), length(gList2), length(intersect(gList1, gList2))))
}