library(Matrix)
library(Seurat)
library(harmony)
library(ggplot2)
library(ggrepel)

## load scTenifoldnet results
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

load('AHR_test/Results/PreenterocytesDR.RData')
O$manifoldAlignment <- O$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O$manifoldAlignment)),]
DR <- dRegulation(O$manifoldAlignment, 'Ahr')
O$diffRegulation <- DR
deZ <- DR$Z
names(deZ) <- toupper(DR$gene)

library(fgsea)
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
deE <- fgsea(REACTOME, deZ, 1e6)

## daniel's method
load('AHR_test/Results/Preenterocytes.RData')
O$manifoldAlignment <- O$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O$manifoldAlignment)),]
head(O$manifoldAlignment)
DR <- dRegulation(O$manifoldAlignment, 'Ahr')
O$diffRegulation <- DR
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)
drE <- fgsea(REACTOME, drZ, 1e6)

library(UpSetR)
png('AHR_test/Results_net/fgsea_daniel_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

## yan's method
O <- readRDS("AHR_test/Results_Yan/Ahr_yan.rds")
O$manifoldAlignment <- O$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O$manifoldAlignment)),]
head(O$manifoldAlignment)
O$manifoldAlignment <- O$manifoldAlignment[, 1:3]
DR <- dRegulation(O$manifoldAlignment, 'Ahr')
O$diffRegulation <- DR
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)
drE <- fgsea(REACTOME, drZ, 1e6)

library(UpSetR)
png('AHR_test/Results_net/fgsea_yan_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

## yan's new method
O <- readRDS("AHR_test/Results_YAN_new/data/Ahr_yan_k2_gamma0_method1.rds")
O$manifoldalignment <- O$manifoldalignment[!grepl('_Rpl|_Rps',rownames(O$manifoldalignment)),]
head(O$manifoldalignment)
O$manifoldalignment <- O$manifoldalignment[, 2:4]
DR <- dRegulation(O$manifoldalignment, 'Ahr')
O$diffRegulation <- DR
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)
drE <- fgsea(REACTOME, drZ, 1e6)

library(UpSetR)
png('AHR_test/Results_net/fgsea_yan_new_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

## yan's new norm method
O <- readRDS('AHR_test/Results_YAN_new_norm/data/Ahr_yan_k2_gamma0_method0.rds')
O$manifoldalignment <- O$manifoldalignment[!grepl('_Rpl|_Rps',rownames(O$manifoldalignment)),]
head(O$manifoldalignment)
O$manifoldalignment <- O$manifoldalignment[, 2:3]
DR <- dRegulation(O$manifoldalignment, 'Ahr')
O$diffRegulation <- DR
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)
drE <- fgsea(REACTOME, drZ, 1e6)

library(UpSetR)
png('AHR_test/Results_net/fgsea_yan_new_norm_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()
