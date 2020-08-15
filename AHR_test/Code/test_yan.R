library(Matrix)
library(scTenifoldNet)

## load the first data set
load('AHR_test/Results/Preenterocytes.RData')

WT <- O$WT
KO <- O$KO

manifoldAlignment <- function(X, Y, d = 5, lambda = 0.9){
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
  W = t(W) %*% W
 
  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
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

dRegulation <- function(manifoldOutput, gKO){
  
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl('^X_', geneList)]
  geneList <- gsub('^X_','', geneList)
  nGenes <- length(geneList)
  
  eGenes <- nrow(manifoldOutput)/2
  
  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl('^Y_', eGeneList)]
  eGeneList <- gsub('^Y_','', eGeneList)
  
  if(nGenes != eGenes){
    stop('Number of identified and expected genes are not the same')
  }
  if(!all(eGeneList == geneList)){
    stop('Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order')
  }
  
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- stats::dist(I)
    O <- as.numeric(O)
    return(O)
  })
  
  ### BOX-COX
  lambdaValues <- seq(-2,2,length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- MASS::boxcox(dMetric~1, plot=FALSE, lambda = lambdaValues)
  BC <- BC$x[which.max(BC$y)]
  if(BC < 0){
    nD <- 1/(dMetric ^ BC)
  } else {
    nD <- dMetric ^ BC
  }
  Z <- scale(nD)
  dOut <- data.frame(
    gene = geneList, 
    distance = dMetric,
    Z = Z
    )
  dOut <- dOut[order(dOut$distance, decreasing = TRUE),]
  FC <- (dOut$distance^2)/mean((dOut$distance[-seq_len(length(gKO))]^2))
  pValues <- pchisq(q = FC,df = 1,lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = 'fdr')
  dOut$FC = FC
  dOut$p.value = pValues
  dOut$p.adj = pAdjusted
  dOut <- as.data.frame.array(dOut)
  return(dOut)
}

set.seed(1)
MA <- manifoldAlignment(WT, KO, d = 5)
set.seed(1)
DR <- dRegulation(MA, gKO = 'Ahr')

out_list <- list()
out_list$WT <- WT
out_list$KO <- KO
out_list$manifoldAlignment <- MA
out_list$diffRegulation <- DR

saveRDS(out_list, "Ahr_yan.rds")

