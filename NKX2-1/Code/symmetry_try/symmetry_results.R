library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')

## inpu data
# load('../Results/GSM3716703.RData')
X_both0 <- readRDS("symmetry_results/MA_both0.rds")
X_both2 <- readRDS("symmetry_results/MA_both2.rds")
X_row0 <- readRDS("symmetry_results/MA_row0.rds")
X_col0 <- readRDS("symmetry_results/MA_col0.rds")

## modify the results
fixPValues <- function(X){
  X <- list("manifoldAlignment" = X)
  X$manifoldAlignment <- X$manifoldAlignment[,1:30]
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-1]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
  return(X)
}

X_both0 <- fixPValues(X_both0)
X_both2 <- fixPValues(X_both2)
X_row0 <- fixPValues(X_row0)
X_col0 <- fixPValues(X_col0)

## marker gene
markerGenes <- read.csv('../Data/pnas.1906663116.sd01.csv', stringsAsFactors = FALSE, row.names = 1)
markerGenes$T.test.p.value <- as.numeric(markerGenes$T.test.p.value)
markerGenes <- markerGenes[complete.cases(markerGenes),]
markerGenesAT1 <- markerGenes$gene_short_name[markerGenes$T.test.p.value < 0.05]

markerGenes <- read.csv('../Data/pnas.1906663116.sd05.csv', stringsAsFactors = FALSE)
markerGenes <- markerGenes[,1:10]
markerGenes$T.test.p.value <- as.numeric(markerGenes$T.test.p.value)
markerGenes <- markerGenes[complete.cases(markerGenes),]
markerGenesAT2 <- markerGenes$gene_short_name[markerGenes$T.test.p.value < 0.05]

markerGenes <- unique(c(markerGenesAT1, markerGenesAT2)) # 7808
length(markerGenes)

######################################################################################################
##################################   row 0 method                   ##################################
######################################################################################################

dGenes <- X_row0$diffRegulation$gene[X_row0$diffRegulation$p.value < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis (failed!)
E <- enrichr(dGenes, c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
E <- do.call(rbind.data.frame, E)
E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))

name_need <- c("GO_Cellular_Component_2018.1", "Reactome_2016.1", "BioPlanet_2019.1",
               "GO_Cellular_Component_2018.4", "GO_Cellular_Component_2018.2")
term_need <- c("Lamellar body", "Surfactant metabolism", "Cell adhesion molecules", "Multivesicular body lumen", "Lysosome")
index1 <- which(rownames(E) %in% name_need)
index2 <- which(E$Term %in% term_need)
index <- intersect(index1, index2)
nrow(E)
index
E <- E[index, c(1,2,3,4)]
E


######################################################################################################
##################################   col 0 method                   ##################################
######################################################################################################

dGenes <- X_col0$diffRegulation$gene[X_col0$diffRegulation$p.value < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis (failed!)
E <- enrichr(dGenes, c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
E <- do.call(rbind.data.frame, E)
E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))

name_need <- c("GO_Cellular_Component_2018.1", "Reactome_2016.1", "BioPlanet_2019.1",
               "GO_Cellular_Component_2018.4", "GO_Cellular_Component_2018.2")
term_need <- c("Lamellar body", "Surfactant metabolism", "Cell adhesion molecules", "Multivesicular body lumen", "Lysosome")
index1 <- which(rownames(E) %in% name_need)
index2 <- which(E$Term %in% term_need)
index <- intersect(index1, index2)
nrow(E)
index
E <- E[index, c(1,2,3,4)]
E

######################################################################################################
##################################   both 0 method                   ##################################
######################################################################################################

dGenes <- X_both0$diffRegulation$gene[X_both0$diffRegulation$p.value < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis (failed!)
E <- enrichr(dGenes, c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
E <- do.call(rbind.data.frame, E)
E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))

name_need <- c("GO_Cellular_Component_2018.1", "Reactome_2016.1", "BioPlanet_2019.1",
               "GO_Cellular_Component_2018.4", "GO_Cellular_Component_2018.2")
term_need <- c("Lamellar body", "Surfactant metabolism", "Cell adhesion molecules", "Multivesicular body lumen", "Lysosome")
index1 <- which(rownames(E) %in% name_need)
index2 <- which(E$Term %in% term_need)
index <- intersect(index1, index2)
nrow(E)
index
E <- E[index, c(1,2,3,4)]
E


######################################################################################################
##################################   both divide by 2  method       ##################################
######################################################################################################

dGenes <- X_both2$diffRegulation$gene[X_both2$diffRegulation$p.value < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis (failed!)
E <- enrichr(dGenes, c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
E <- do.call(rbind.data.frame, E)
E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))

name_need <- c("GO_Cellular_Component_2018.1", "Reactome_2016.1", "BioPlanet_2019.1",
               "GO_Cellular_Component_2018.4", "GO_Cellular_Component_2018.2")
term_need <- c("Lamellar body", "Surfactant metabolism", "Cell adhesion molecules", "Multivesicular body lumen", "Lysosome")
index1 <- which(rownames(E) %in% name_need)
index2 <- which(E$Term %in% term_need)
index <- intersect(index1, index2)
E <- E[index, c(1,2,3,4)]
E
