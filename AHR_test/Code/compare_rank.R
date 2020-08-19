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
load('AHR_test/Results/PreenterocytesDR.RData')
O_net <- O
O_net$manifoldAlignment <- O_net$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O_net$manifoldAlignment)),]
DR <- dRegulation(O_net$manifoldAlignment, 'Ahr')
O_net$diffRegulation <- DR
Dr_net <- O_net$diffRegulation[, 1:2]
Dr_net$order <- 1:nrow(Dr_net)
rownames(Dr_net) <- Dr_net$gene
Dr_net <- Dr_net[, 2:3]

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

## yan's method
O_yan <- readRDS("AHR_test/Results_Yan/Ahr_yan.rds")
O_yan$manifoldAlignment <- O_yan$manifoldAlignment[!grepl('_Rpl|_Rps',rownames(O_yan$manifoldAlignment)),]
head(O_yan$manifoldAlignment)
O_yan$manifoldAlignment <- O_yan$manifoldAlignment[, 1:3]
DR <- dRegulation(O_yan$manifoldAlignment, 'Ahr')
O_yan$diffRegulation <- DR
Dr_yan <- O_yan$diffRegulation[, 1:2]
Dr_yan$order <- 1:nrow(Dr_yan)
rownames(Dr_yan) <- Dr_yan$gene
Dr_yan <- Dr_yan[, 2:3]

load("AHR_test/Results/DR_ori.Rdata")
Dr_de <- DE[, 2]
Dr_de <- as.data.frame(Dr_de)
Dr_de$gene <- rownames(DE)
Dr_de$order <- 1:nrow(DE)
colnames(Dr_de) <- c("FC", "gene", "order")
rownames(Dr_de) <- Dr_de$gene
Dr_de <- Dr_de[, c(1, 3)]

########### plot ###############

## yan and daniel
gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_yan)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_yan_tmp <- Dr_yan[g_shared, ]
png('AHR_test/Compare_gene_list/order_daniel_yan.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = Dr_daniel_tmp$order, y = Dr_yan_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "daniel order", y = "yan order", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_yan)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_yan_tmp <- Dr_yan[g_shared, ]
png('AHR_test/Compare_gene_list/distance_daniel_yan.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = log(Dr_daniel_tmp$distance), y = log(Dr_yan_tmp$distance))) +
  geom_point(size = 0.1) +
  labs(x = "daniel log distance", y = "yan log distance", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

## yan and net
gList1 <- rownames(Dr_net)
gList2 <- rownames(Dr_yan)
g_shared <- intersect(gList1, gList2)
Dr_net_tmp <- Dr_net[g_shared, ]
Dr_yan_tmp <- Dr_yan[g_shared, ]
png('AHR_test/Compare_gene_list/order_net_yan.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = Dr_net_tmp$order, y = Dr_yan_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "net order", y = "yan order", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

gList1 <- rownames(Dr_net)
gList2 <- rownames(Dr_yan)
g_shared <- intersect(gList1, gList2)
Dr_net_tmp <- Dr_net[g_shared, ]
Dr_yan_tmp <- Dr_yan[g_shared, ]
png('AHR_test/Compare_gene_list/distance_net_yan.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = log(Dr_net_tmp$distance), y = log(Dr_yan_tmp$distance))) +
  geom_point(size = 0.1) +
  labs(x = "net log distance", y = "yan log distance", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

## daniel and net
gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_net)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_net_tmp <- Dr_net[g_shared, ]
png('AHR_test/Compare_gene_list/order_daniel_net.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = Dr_daniel_tmp$order, y = Dr_net_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "daniel order", y = "net order", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

gList1 <- rownames(Dr_daniel)
gList2 <- rownames(Dr_net)
g_shared <- intersect(gList1, gList2)
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
Dr_net_tmp <- Dr_net[g_shared, ]
png('AHR_test/Compare_gene_list/distance_daniel_net.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = log(Dr_daniel_tmp$distance), y = log(Dr_net_tmp$distance))) +
  geom_point(size = 0.1) +
  labs(x = "daniel log distance", y = "net log distance", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

######################## Compare with DE ########################
## yan and de
gList1 <- rownames(Dr_de)
gList2 <- rownames(Dr_yan)
g_shared <- intersect(gList1, gList2)
Dr_de_tmp <- Dr_de[g_shared, ]
Dr_yan_tmp <- Dr_yan[g_shared, ]
png('AHR_test/Compare_gene_list/order_de_yan.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = Dr_de_tmp$order, y = Dr_yan_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "de order", y = "yan order", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

## daniel and de
gList1 <- rownames(Dr_de)
gList2 <- rownames(Dr_daniel)
g_shared <- intersect(gList1, gList2)
Dr_de_tmp <- Dr_de[g_shared, ]
Dr_daniel_tmp <- Dr_daniel[g_shared, ]
png('AHR_test/Compare_gene_list/order_de_daniel.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = Dr_de_tmp$order, y = Dr_daniel_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "de order", y = "daniel order", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

## net and de
gList1 <- rownames(Dr_de)
gList2 <- rownames(Dr_net)
g_shared <- intersect(gList1, gList2)
Dr_de_tmp <- Dr_de[g_shared, ]
Dr_net_tmp <- Dr_net[g_shared, ]
png('AHR_test/Compare_gene_list/order_de_net.png', width = 1200, height = 1200, res = 300)
ggplot(data = NULL, aes(x = Dr_de_tmp$order, y = Dr_net_tmp$order)) +
  geom_point(size = 0.1) +
  labs(x = "de order", y = "net order", title = paste0("Total ", length(g_shared), " gene"))
dev.off()

######################## Compare with DE ########################
## yan's method
png('AHR_test/Compare_gene_list/order_distance_yan.png', width = 1200, height = 1200, res = 300)
ggplot(data = Dr_yan, aes(x = order, y = log(distance))) +
  geom_point(size = 0.1) +
  labs(x = "order", y = "log distance", title = paste0("Yan: total ", length(g_shared), " gene"))
dev.off()

## daniel's method
png('AHR_test/Compare_gene_list/order_distance_daniel.png', width = 1200, height = 1200, res = 300)
ggplot(data = Dr_daniel, aes(x = order, y = log(distance))) +
  geom_point(size = 0.1) +
  labs(x = "order", y = "log distance", title = paste0("Daniel: total ", length(g_shared), " gene"))
dev.off()

## net method
png('AHR_test/Compare_gene_list/order_distance_net.png', width = 1200, height = 1200, res = 300)
ggplot(data = Dr_net, aes(x = order, y = log(distance))) +
  geom_point(size = 0.1) +
  labs(x = "order", y = "log distance", title = paste0("Net: total ", length(g_shared), " gene"))
dev.off()




