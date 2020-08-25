library(Matrix)
library(Seurat)
library(fgsea)
library(UpSetR)
library(OrderedList)

mmuKEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')

#### DE ####
WT <- Read10X_h5('MALAT1_test/data/WT.h5')
WT <- as.matrix(WT)
genesWT <- readLines('MALAT1_test/data/genesWT.txt')
barcodesWT <- readLines('MALAT1_test/data/barcodesWT.txt')
WT <- WT[genesWT, ]
WT <- WT[, barcodesWT]
WT <-WT[!rownames(WT) %in% 'Malat1',]

KO <- Read10X_h5('MALAT1_test/data/KO.h5')
genesKO <- readLines('MALAT1_test/data/genesKO.txt')
barcodesKO <- readLines('MALAT1_test/data/barcodesKO.txt')
KO <- KO[genesKO, ]
KO <- KO[, barcodesKO]
KO <- KO[!rownames(KO) %in% 'Malat1',]

WT <- CreateSeuratObject(WT, project = 'WT')
KO <- CreateSeuratObject(KO, project = 'KO')

ALL <- merge(WT,KO)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL, verbose = FALSE)
ALL <- RunPCA(ALL, verbose = FALSE)
ALL <- RunUMAP(ALL, dims = 1:50, verbose = FALSE)
UMAPPlot(ALL)
DE <- FindMarkers(ALL, ident.1 = 'KO', ident.2 = 'WT', test.use = 'MAST', logfc.threshold = 0)
FC <- DE$avg_logFC
names(FC) <- toupper(rownames(DE))

#### DR ####
load('betaMALATko.RData')
DR <- MALAT1$diffRegulation
Z <- DR$Z
names(Z) <- toupper(DR$gene)


#### Comparison ####
sGenes <- intersect(names(Z), names(FC))
Z <- Z[sGenes]
FC <- abs(FC[sGenes])

eDR <- fgseaMultilevel(mmuKEGG, Z)
eDE <- fgseaMultilevel(mmuKEGG, FC)

eDR <- eDR[eDR$NES > 0 & eDR$padj < 0.05]
eDE <- eDE[eDE$NES > 0 & eDE$padj < 0.05]

upset(fromList(list(DE=eDE$pathway, DR=eDR$pathway)))

lZ <- names(sort(Z, decreasing = TRUE))
lFC <- names(sort(FC, decreasing = TRUE))

plot(compareLists(lZ,lFC, no.reverse = TRUE, alphas = c(0.005)))

