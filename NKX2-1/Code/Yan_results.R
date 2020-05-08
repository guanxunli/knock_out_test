library(fgsea)
library(ggplot2)
library(enrichR)
library(igraph)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')

## inpu data
X1_row0 <- readRDS("Yan_results/Yan_row0.rds")
X1_col0 <- readRDS("Yan_results/Yan_col0.rds")
X1_both0 <- readRDS("Yan_results/Yan_both0.rds")
X1_both2 <- readRDS("Yan_results/Yan_both2.rds")

X2_row0 <- readRDS("Yan_results/Yan2_row0.rds")
X2_col0 <- readRDS("Yan_results/Yan2_col0.rds")
X2_both0 <- readRDS("Yan_results/Yan2_both0.rds")
X2_both2 <- readRDS("Yan_results/Yan2_both2.rds")

## modify the results
fixPValues <- function(X){
  X <- list("manifoldAlignment" = X)
  X$manifoldAlignment <- X$manifoldAlignment[,1:50]
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-1]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
  return(X)
}

## make correct
X1_row0 <- fixPValues(X1_row0)
X1_col0 <- fixPValues(X1_col0)
X1_both0 <- fixPValues(X1_both0)
X1_both2 <- fixPValues(X1_both2)
X2_row0 <- fixPValues(X2_row0)
X2_col0 <- fixPValues(X2_col0)
X2_both0 <- fixPValues(X2_both0)
X2_both2 <- fixPValues(X2_both2)


# marker gene
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
################################## method1 row 0                #####################################
######################################################################################################

## get gene list
dGenes <- X1_row0$diffRegulation$gene[X1_row0$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis
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
################################## method1 col 0                #####################################
######################################################################################################

## get gene list
dGenes <- X1_col0$diffRegulation$gene[X1_col0$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis
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
################################## method1 both 0                #####################################
######################################################################################################

## get gene list
dGenes <- X1_both0$diffRegulation$gene[X1_both0$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis
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

######################################################################################################
################################## method1 both 2                #####################################
######################################################################################################

## get gene list
dGenes <- X1_both2$diffRegulation$gene[X1_both2$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis
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
##################################      method2 row 0                #####################################
######################################################################################################

## get gene list
dGenes <- X2_row0$diffRegulation$gene[X2_row0$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis
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
##################################         method2 col 0             #################################
######################################################################################################

## get gene list
dGenes <- X2_col0$diffRegulation$gene[X2_col0$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))


## enrichment analysis
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
##################################     method2 both 0            #####################################
######################################################################################################

## get gene list
dGenes <- X2_both0$diffRegulation$gene[X2_both0$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis
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
##################################     method2 both 2            #####################################
######################################################################################################

## get gene list
dGenes <- X2_both2$diffRegulation$gene[X2_both2$diffRegulation$p.adj < 0.05]
length(dGenes)
length(intersect(markerGenes, dGenes))

## enrichment analysis
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

