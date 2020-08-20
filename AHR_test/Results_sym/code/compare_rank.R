library(Matrix)
library(Seurat)
library(harmony)
library(ggplot2)
library(ggrepel)

## Add function
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

####################################### Gene compare ranke ##############################################
## daniel's method
load('AHR_test/Results/Preenterocytes.RData')
O_daniel <- O 
O_daniel$manifoldAlignment <- O_daniel$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O_daniel$manifoldAlignment)),]
head(O_daniel$manifoldAlignment)
DR <- dRegulation(O_daniel$manifoldAlignment, 'Ahr')
O_daniel$diffRegulation <- DR
Dr_daniel <- O_daniel$diffRegulation[, 1:2]
Dr_daniel$order <- 1:nrow(Dr_daniel)
rownames(Dr_daniel) <- Dr_daniel$gene
Dr_daniel <- Dr_daniel[, 2:3]

## daniel1
MA <- readRDS("AHR_test/Results_sym/data/ahr_daniel1.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
Dr_daniel1 <- data.frame("distance" = DR$distance, "order" = 1:nrow(DR))
rownames(Dr_daniel1) <- DR$gene

## daniel2
MA <- readRDS("AHR_test/Results_sym/data/ahr_daniel2.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
Dr_daniel2 <- data.frame("distance" = DR$distance, "order" = 1:nrow(DR))
rownames(Dr_daniel2) <- DR$gene

## daniel3
MA <- readRDS("AHR_test/Results_sym/data/ahr_daniel3.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
Dr_daniel3 <- data.frame("distance" = DR$distance, "order" = 1:nrow(DR))
rownames(Dr_daniel3) <- DR$gene

## yan1
MA <- readRDS("AHR_test/Results_sym/data/ahr_yan1.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
Dr_yan1 <- data.frame("distance" = DR$distance, "order" = 1:nrow(DR))
rownames(Dr_yan1) <- DR$gene

########### plot ###############

## daniel and daniel1
gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_daniel1)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_daniel1_tmp <- Dr_daniel1[g_shared, ]
ggplot(data = NULL, aes(x = Dr_daniel_tmp$order, y = Dr_daniel1_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "daniel order", y = "daniel1 order", title = paste0("Total ", length(g_shared), " gene"))

## daniel and daniel2
gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_daniel2)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_daniel2_tmp <- Dr_daniel2[g_shared, ]
ggplot(data = NULL, aes(x = Dr_daniel_tmp$order, y = Dr_daniel2_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "daniel order", y = "daniel2 order", title = paste0("Total ", length(g_shared), " gene"))

## daniel and daniel3
gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_daniel3)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_daniel3_tmp <- Dr_daniel3[g_shared, ]
ggplot(data = NULL, aes(x = Dr_daniel_tmp$order, y = Dr_daniel3_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "daniel order", y = "daniel3 order", title = paste0("Total ", length(g_shared), " gene"))

## daniel and yan1
gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_yan1)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_yan1_tmp <- Dr_yan1[g_shared, ]
ggplot(data = NULL, aes(x = Dr_daniel_tmp$order, y = Dr_yan1_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "daniel order", y = "yan1 order", title = paste0("Total ", length(g_shared), " gene"))

## daniel11 and daniel2
gList1 <- rownames(Dr_daniel1)
gList2 <- rownames(Dr_daniel2)
g_shared <- intersect(gList1, gList2)
Dr_daniel1_tmp <- Dr_daniel1[g_shared, ]
Dr_daniel2_tmp <- Dr_daniel2[g_shared, ]
ggplot(data = NULL, aes(x = Dr_daniel1_tmp$order, y = Dr_daniel2_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "daniel1 order", y = "daniel2 order", title = paste0("Total ", length(g_shared), " gene"))


