library(Matrix)
library(Seurat)
library(harmony)
library(ggplot2)
library(ggrepel)

load("AHR_test/Results/DR_ori.Rdata")

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

## daniel1
MA <- readRDS("AHR_test/Results_sym/data/ahr_daniel1.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

# drZ <- DR$distance
# names(drZ) <- toupper(DR$gene)

library(fgsea)
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
drE <- fgsea(REACTOME, drZ, 1e6)
deE <- fgsea(REACTOME, deZ, 1e6)

library(UpSetR)
png('AHR_test/Results_sym/results/fgsea_daniel1.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

## daniel2
MA <- readRDS("AHR_test/Results_sym/data/ahr_daniel2.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

# drZ <- DR$distance
# names(drZ) <- toupper(DR$gene)

library(fgsea)
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
drE <- fgsea(REACTOME, drZ, 1e6)

library(UpSetR)
png('AHR_test/Results_sym/results/fgsea_daniel2.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

## daniel3
MA <- readRDS("AHR_test/Results_sym/data/ahr_daniel3.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

# drZ <- DR$distance
# names(drZ) <- toupper(DR$gene)

library(fgsea)
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
drE <- fgsea(REACTOME, drZ, 1e6)

library(UpSetR)
png('AHR_test/Results_sym/results/fgsea_daniel3.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

## yan1
MA <- readRDS("AHR_test/Results_sym/data/ahr_yan1.rds")
head(MA)
MA <- MA[, 1:2]
DR <- dRegulation(MA, 'Ahr')
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

# drZ <- DR$distance
# names(drZ) <- toupper(DR$gene)

library(fgsea)
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
drE <- fgsea(REACTOME, drZ, 1e6)

library(UpSetR)
png('AHR_test/Results_sym/results/fgsea_yan1.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()
