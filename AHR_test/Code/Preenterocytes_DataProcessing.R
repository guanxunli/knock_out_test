library(Matrix)
library(Seurat)
library(harmony)
library(ggplot2)
library(ggrepel)

load('AHR_test/Results/Preenterocytes.RData')
O$manifoldAlignment <- O$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O$manifoldAlignment)),]
DR <- scTenifoldKnk:::dRegulation(O$manifoldAlignment, 'Ahr')
O$diffRegulation <- DR

WT <- readMM('AHR_test/Data/Preenterocytes_WT.mtx')
rownames(WT) <- readLines('AHR_test/Data/Preenterocytes_WT_genes.txt')
colnames(WT) <- readLines('AHR_test/Data/Preenterocytes_WT_barcodes.txt')
WT <- WT[!grepl('^MT-|^RPL|^RPS',rownames(WT), ignore.case = TRUE),]

  KO <- readMM('AHR_test/Data/Preenterocytes_KO.mtx')
rownames(KO) <- readLines('AHR_test/Data/Preenterocytes_KO_genes.txt')
colnames(KO) <- readLines('AHR_test/Data/Preenterocytes_KO_barcodes.txt')
KO <- KO[!grepl('^MT-|^RPL|^RPS',rownames(KO), ignore.case = TRUE),]

WT <- CreateSeuratObject(WT, project = 'WT')
KO <- CreateSeuratObject(KO, project = 'KO')

ALL <- merge(WT,KO)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- RunPCA(ALL)
ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
ALL <- RunTSNE(ALL, reduction = 'harmony')
ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:50)
ALL <- RunUMAP(ALL, reduction = 'pca', dims = 1:50)
UMAPPlot(ALL)

DE <- FindMarkers(ALL, ident.1 = 'KO', ident.2 = 'WT', test.use = 'MAST', logfc.threshold = 0)
DE <- DE[order(abs(DE$avg_logFC), decreasing = TRUE),]
deList <- rownames(DE)
deZ <- abs(DE$avg_logFC)
names(deZ) <- toupper(rownames(DE))

dF <- data.frame(FC=DE$avg_logFC, P=-log10(DE$p_val), G = rownames(DE))
dF <- as.data.frame.array(dF)
dF$O <- seq_len(nrow(dF))
dF <- dF[order(abs(dF$FC)),]
dF$O[!is.na(dF$G)] <- 1e6
dF <- dF[order(dF$O, decreasing = FALSE),]
dF$COL <- densCols(dF[,1:2])
dF$COL[dF$G %in% DR$gene[DR$p.adj < 0.05]] <- 'red'
dF$G[!dF$G %in% DR$gene[DR$p.adj < 0.05]] <- NA
gT <- ifelse(is.na(dF$G),0.2,1)
BG <- rownames(DE)[DE$p_val_adj < 0.05 & abs(DE$avg_logFC) > 0.1]
BG <- intersect(dF$G, BG)
FF <- ifelse(dF$G %in% BG, 2, 1)

png('VLN3.png', width = 1500, height = 1500, res=300)
ggplot(dF, aes(FC, P, label = G)) + geom_point(color = dF$COL, alpha = gT) + xlim(c(-0.55,0.55)) + theme_bw() + xlab(log[2]~'Fold-Change') + ylab(-log[10]~'P-value') + geom_text_repel(box.padding = 0.15, segment.size = 0.05, aes(fontface= FF))
dev.off()


drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

drList <- DR$gene
deList <- deList[deList %in% drList]
drList <- drList[drList %in% deList]

# CLA <- compareList(deList,drList,B = 100, label = 'MAST (abs FC) vs\nscTenifoldKnk')
# png('cl1_Ahr.png', width = 1200, height = 1200, res = 300)
# CLA
# dev.off()

# library(OrderedList)
# png('DE_DR2.png', width = 2000, height = 1000, res = 300)
# par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
# plot(compareLists(deList,drList, alphas = 0.005))
# dev.off()


library(fgsea)
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
drE <- fgsea(REACTOME, drZ, 1e6)
deE <- fgsea(REACTOME, deZ, 1e6)

library(UpSetR)
png('fgsea3.png', width = 600, height = 600, res = 300)
upset(fromList(list(Simulation=drE$pathway[drE$padj < 0.05 & drE$NES > 0],Real=deE$pathway[deE$padj < 0.05 & deE$NES > 0])))
dev.off()

library(ggplot2)
library(ggrepel)
pData <- data.frame(FC=DE$avg_logFC, P=-log10(DE$p_val), G=rownames(DE))
pData$G[!pData$G %in% DR$gene[DR$p.adj < 0.05]] <- NA
pData$L <- seq_len(nrow(pData))
pData$L[!is.na(pData$G)] <- 1e6
pData <- pData[order(pData$L, decreasing = FALSE),]
pData$C <- densCols(pData[,1:2])
pData$C[!is.na(pData$G)] <- 'red'
png('DE_drLabels.png', width = 2000, height = 2000, res = 300)
ggplot(pData, aes(FC,P,label=G)) + geom_point(color=pData$C) + geom_text_repel() + theme_bw() + ylim(c(0,80)) + xlim(c(-.6,.6))
dev.off()       

source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')

png('KO3.png', width = 3500, height = 3500, res = 300)
plotKO(O, 'Ahr', nCategories = 10)
dev.off()




#### KNK vs NET #### 
load('../Results/Preenterocytes.RData')
MA <- O$manifoldAlignment
MA <- MA[!grepl('_MT-|_RPL|_RPS',rownames(MA), ignore.case = TRUE),]
DR <- scTenifoldNet::dRegulation(MA)
#DR <- O$diffRegulation
DR$FC <- (DR$distance^2)/mean(DR$distance[-1]^2)
DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
O$diffRegulation <- DR
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

load('../Results/PreenterocytesDR.RData')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
png('KNK_NET.png', width = 1500, height = 1500, res = 300)
plotDR(O, labelGenes = intersect(DR$gene[DR$p.adj < 0.05], O$diffRegulation$gene[O$diffRegulation$p.value < 0.05]), boldGenes = intersect(DR$gene[DR$p.adj < 0.05], O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]))
dev.off()

realKO <- O$diffRegulation
realZ <- realKO$FC
names(realZ) <- toupper(realKO$gene)
realList <- realKO$gene
drList <- DR$gene
realList <- realList[realList %in% drList]
#drList <- drList[drList %in% realList]

library(OrderedList)
library(ggplot2)
library(ggrepel)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/compareLists.R')
CLB <- compareList(realList,drList,B = 100, label = 'scTenifoldNet vs\nscTenifoldKnk')
png('cl2_Ahr.png', width = 1200, height = 1200, res = 300)
CLB
dev.off()
#plot(compareLists(realList,drList, alphas = 0.005))

rownames(DR) <- DR$gene
rownames(realKO) <- realKO$gene

DR <- DR[realList,]
realKO <- realKO[realList,]

dF <- data.frame(scTenifoldNet=realKO$Z, scTenifoldKnk = DR$Z)
dF$G <- realList
dF$G[DR$p.value > 0.05] <- NA
dF$G[realKO$p.value > 0.05] <- NA
gColor <- densCols(dF[,1:2])
gColor[!is.na(dF$G)] <- 'red'
png('KNK_NET.png', height = 2000, width = 2000, res = 300)
ggplot(dF, aes(scTenifoldNet, scTenifoldKnk, label= G)) + theme_bw() + geom_point(col= gColor) + geom_text_repel()
dev.off()

# library(fgsea)
# KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
# REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
# sGenes <- intersect(names(drZ), names(realZ))
# drE <- fgsea(REACTOME, drZ[sGenes], 1e3)
# koE <- fgsea(REACTOME, realZ[sGenes], 1e3)
# library(UpSetR)
# #png('fgsea.png', width = 600, height = 600, res = 300)
# upset(fromList(list(Simulation=drE$pathway[drE$pval < 0.05 & drE$NES > 0],Real=koE$pathway[koE$pval < 0.05 & koE$NES > 0])))
# #dev.off()


library(fgsea)
load('../Results/Preenterocytes.RData')
MA <- O$manifoldAlignment
MA <- MA[!grepl('_MT-|_RPL|_RPS',rownames(MA), ignore.case = TRUE),]
DR <- scTenifoldNet::dRegulation(MA)
#DR <- O$diffRegulation
DR$FC <- (DR$distance^2)/mean(DR$distance[-1]^2)
DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
O$diffRegulation <- DR
drZ <- DR$Z
names(drZ) <- toupper(DR$gene)

# drZ <- DR$Z
# names(drZ) <- toupper(DR$gene)
# CT <- read.table('PanglaoDB_markers_27_Mar_2020.tsv.gz', sep = '\t', header = TRUE)
# CT <- CT[CT$organ %in% 'GI tract',]
CT <- clustermole::clustermole_markers()
CT <- CT[CT$species %in% 'Mouse',]
CT <- CT[grepl('Intest', CT$organ, ignore.case = TRUE),]
ctNames <- unique(CT$celltype)
CT <- lapply(ctNames, function(X){
  unique(CT$gene[CT$celltype %in% X])
})
names(CT) <- ctNames

# CM <- gmtPathways('../Data/CannonicalMarkers.txt')
# CM <- lapply(CM, toupper)
library(ggplot2)
png('DE_gseaMarkers.png', width = 1000, height = 1000, res = 300)
E <- fgseaMultilevel(CT, deZ)
plotEnrichment(CT$Enterocyte, deZ) + xlab('Gene rank') + ylab('Enrichment Score') + labs(title = 'Enterocytes', subtitle = paste0('MAST\nFDR = ', formatC(E$padj[E$pathway == 'Enterocyte'], 2, format = 'e'))) + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png('DR_gseaMarkers.png', width = 1000, height = 1000, res = 300)
E <- fgseaMultilevel(CT, drZ)
plotEnrichment(CT$Enterocyte, drZ) + xlab('Gene rank') + ylab('Enrichment Score') + labs(title = 'Enterocytes', subtitle = paste0('scTenifoldKnk\nFDR = ', formatC(E$padj[E$pathway == 'Enterocyte'], 2, format = 'e'))) + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()