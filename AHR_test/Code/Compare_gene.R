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

# load("AHR_test/Results/DR_ori.Rdata")
# gene_DE <- rownames(DE)[DE$p_val_adj < 0.05]
# length(gene_DE)

####################################### Gene compare ##############################################
load('AHR_test/Results/PreenterocytesDR.RData')
O$manifoldAlignment <- O$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O$manifoldAlignment)),]
DR <- dRegulation(O$manifoldAlignment, 'Ahr')
O$diffRegulation <- DR
gList_net <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]

## daniel's method
load('AHR_test/Results/Preenterocytes.RData')
O$manifoldAlignment <- O$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O$manifoldAlignment)),]
head(O$manifoldAlignment)
DR <- dRegulation(O$manifoldAlignment, 'Ahr')
O$diffRegulation <- DR
gList_daniel <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]

library(UpSetR)
png('AHR_test/Compare_gene_list/glist_daniel_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(scTenifoldNet=gList_net, daniel_knk=gList_daniel)))
dev.off()

## yan's method
O <- readRDS("AHR_test/Results_Yan/Ahr_yan.rds")
O$manifoldAlignment <- O$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O$manifoldAlignment)),]
head(O$manifoldAlignment)
O$manifoldAlignment <- O$manifoldAlignment[, 1:3]
DR <- dRegulation(O$manifoldAlignment, 'Ahr')
O$diffRegulation <- DR
gList_yan <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]

png('AHR_test/Compare_gene_list/glist_yan_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(scTenifoldNet=gList_net, yan_knk=gList_yan)))
dev.off()

png('AHR_test/Compare_gene_list/glist_yan_daniel.png', width = 600, height = 600, res = 300)
upset(fromList(list(daniel_knk=gList_daniel, yan_knk=gList_yan)))
dev.off()

## yan's new method
O <- readRDS("AHR_test/Results_YAN_new/data/Ahr_yan_k2_gamma0_method1.rds")
O$manifoldalignment <- O$manifoldalignment[!grepl('_Rpl|_Rps',rownames(O$manifoldalignment)),]
head(O$manifoldalignment)
O$manifoldalignment <- O$manifoldalignment[, 2:4]
DR <- dRegulation(O$manifoldalignment, 'Ahr')
O$diffRegulation <- DR
gList_yan_new <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]

png('AHR_test/Compare_gene_list/glist_yan_net_new.png', width = 600, height = 600, res = 300)
upset(fromList(list(scTenifoldNet=gList_net, yan_knk=gList_yan_new)))
dev.off()

## yan's new norm method
O <- readRDS('AHR_test/Results_YAN_new_norm/data/Ahr_yan_k2_gamma0_method0.rds')
O$manifoldalignment <- O$manifoldalignment[!grepl('_Rpl|_Rps',rownames(O$manifoldalignment)),]
head(O$manifoldalignment)
O$manifoldalignment <- O$manifoldalignment[, 2:3]
DR <- dRegulation(O$manifoldalignment, 'Ahr')
gList_yan_new_norm <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]

png('AHR_test/Compare_gene_list/glist_yan_net_new_norm.png', width = 600, height = 600, res = 300)
upset(fromList(list(scTenifoldNet=gList_net, yan_knk=gList_yan_new_norm)))
dev.off()

############################## significiant gene pathway ##############################################
library(enrichR)
MET_net <- enrichr(gList_net, c('Reactome_2016','BioPlanet_2019','KEGG_2019_Mouse', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018'))
MET_net <- do.call(rbind.data.frame, MET_net)
MET_net <- MET_net[order(MET_net$Adjusted.P.value),]
MET_net <- MET_net[MET_net$Adjusted.P.value < 0.05,]
pathway_net <- unique(MET_net$Term)

MET_daniel <- enrichr(gList_daniel, c('Reactome_2016','BioPlanet_2019','KEGG_2019_Mouse', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018'))
MET_daniel <- do.call(rbind.data.frame, MET_daniel)
MET_daniel <- MET_daniel[order(MET_daniel$Adjusted.P.value),]
MET_daniel <- MET_daniel[MET_daniel$Adjusted.P.value < 0.05,]
pathway_daniel <- unique(MET_daniel$Term)

MET_yan <- enrichr(gList_yan, c('Reactome_2016','BioPlanet_2019','KEGG_2019_Mouse', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018'))
MET_yan <- do.call(rbind.data.frame, MET_yan)
MET_yan <- MET_yan[order(MET_yan$Adjusted.P.value),]
MET_yan <- MET_yan[MET_yan$Adjusted.P.value < 0.05,]
pathway_yan <- unique(MET_yan$Term)

MET_yan_new <- enrichr(gList_yan_new, c('Reactome_2016','BioPlanet_2019','KEGG_2019_Mouse', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018'))
MET_yan_new <- do.call(rbind.data.frame, MET_yan_new)
MET_yan_new <- MET_yan_new[order(MET_yan_new$Adjusted.P.value),]
MET_yan_new <- MET_yan_new[MET_yan_new$Adjusted.P.value < 0.05,]
pathway_yan_new <- unique(MET_yan_new$Term)

MET_yan_new_norm <- enrichr(gList_yan_new_norm, c('Reactome_2016','BioPlanet_2019','KEGG_2019_Mouse', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018'))
MET_yan_new_norm <- do.call(rbind.data.frame, MET_yan_new_norm)
MET_yan_new_norm <- MET_yan_new_norm[order(MET_yan_new_norm$Adjusted.P.value),]
MET_yan_new_norm <- MET_yan[MET_yan_new_norm$Adjusted.P.value < 0.05,]
pathway_yan_new_norm <- unique(MET_yan_new_norm$Term)

library(UpSetR)
png('AHR_test/Compare_gene_list/pathway_daniel_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(scTenifoldNet=pathway_net, daniel_knk=pathway_daniel)))
dev.off()

png('AHR_test/Compare_gene_list/pathway_yan_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(scTenifoldNet=pathway_net, yan_knk=pathway_yan)))
dev.off()

png('AHR_test/Compare_gene_list/pathway_daniel_yan.png', width = 600, height = 600, res = 300)
upset(fromList(list(daniel_knk=pathway_daniel, yan_knk=pathway_yan)))
dev.off()

# png('AHR_test/Compare_gene_list/pathway_yan_new_net.png', width = 600, height = 600, res = 300)
# upset(fromList(list(scTenifoldNet=pathway_net, yan_knk=pathway_yan_new)))
# dev.off()

png('AHR_test/Compare_gene_list/pathway_yan_new_norm_net.png', width = 600, height = 600, res = 300)
upset(fromList(list(scTenifoldNet=pathway_net, daniel_knk=pathway_yan_new_norm)))
dev.off()

