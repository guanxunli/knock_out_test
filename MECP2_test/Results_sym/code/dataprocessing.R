library(UpSetR)
library(pbapply)
library(ggplot2)
library(enrichR)
library(igraph)
library(fgsea)
library(statsExpressions)
library(patchwork)
library(UpSetR)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')

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

## daniel
SRS3059998 <- readRDS("MECP2_test/Results_Daniel/SRS3059998_daniel.rds")
SRS3059999 <- readRDS("MECP2_test/Results_Daniel/SRS3059999_daniel.rds")

gList1 <- SRS3059998$diffRegulation$gene[SRS3059998$diffRegulation$p.adj < 0.05]
gList2 <- SRS3059999$diffRegulation$gene[SRS3059999$diffRegulation$p.adj < 0.05]

upset(fromList(list(SRS3059998 = gList1, SRS3059999 = gList2)), text.scale = 1.5)

## daniel1
MA1 <- readRDS("MECP2_test/Results_sym/data/daniel1_SRS3059998.rds")
MA2 <- readRDS("MECP2_test/Results_sym/data/daniel1_SRS3059999.rds")
MA1 <- MA1[, 1:2]
MA2 <- MA2[, 1:2]
DR1 <- dRegulation(MA1, "Mecp2")
DR2 <- dRegulation(MA2, "Mecp2")
gList1 <- DR1$gene[DR1$p.adj < 0.05]
gList2 <- DR2$gene[DR2$p.adj < 0.05]

upset(fromList(list(SRS3059998 = gList1, SRS3059999 = gList2)), text.scale = 1.5)

## daniel2
MA1 <- readRDS("MECP2_test/Results_sym/data/daniel2_SRS3059998.rds")
MA2 <- readRDS("MECP2_test/Results_sym/data/daniel2_SRS3059999.rds")
MA1 <- MA1[, 1:2]
MA2 <- MA2[, 1:2]
DR1 <- dRegulation(MA1, "Mecp2")
DR2 <- dRegulation(MA2, "Mecp2")
gList1 <- DR1$gene[DR1$p.adj < 0.05]
gList2 <- DR2$gene[DR2$p.adj < 0.05]

upset(fromList(list(SRS3059998 = gList1, SRS3059999 = gList2)), text.scale = 1.5)

## daniel3
MA1 <- readRDS("MECP2_test/Results_sym/data/daniel3_SRS3059998.rds")
MA2 <- readRDS("MECP2_test/Results_sym/data/daniel3_SRS3059999.rds")
MA1 <- MA1[, 1:2]
MA2 <- MA2[, 1:2]
DR1 <- dRegulation(MA1, "Mecp2")
DR2 <- dRegulation(MA2, "Mecp2")
gList1 <- DR1$gene[DR1$p.adj < 0.05]
gList2 <- DR2$gene[DR2$p.adj < 0.05]

upset(fromList(list(SRS3059998 = gList1, SRS3059999 = gList2)), text.scale = 1.5)

## yan
SRS3059998 <- readRDS("MECP2_test/Results_yan/yan_SRS3059998.rds")
SRS3059999 <- readRDS("MECP2_test/Results_yan/yan_SRS3059999.rds")

gList1 <- SRS3059998$diffRegulation$gene[SRS3059998$diffRegulation$p.adj < 0.05]
gList2 <- SRS3059999$diffRegulation$gene[SRS3059999$diffRegulation$p.adj < 0.05]

upset(fromList(list(SRS3059998 = gList1, SRS3059999 = gList2)), text.scale = 1.5)

## yan1
MA1 <- readRDS("MECP2_test/Results_sym/data/yan1_SRS3059998.rds")
MA2 <- readRDS("MECP2_test/Results_sym/data/yan1_SRS3059999.rds")
DR1 <- dRegulation(MA1, "Mecp2")
DR2 <- dRegulation(MA2, "Mecp2")
gList1 <- DR1$gene[DR1$p.adj < 0.05]
gList2 <- DR2$gene[DR2$p.adj < 0.05]

upset(fromList(list(SRS3059998 = gList1, SRS3059999 = gList2)), text.scale = 1.5)

