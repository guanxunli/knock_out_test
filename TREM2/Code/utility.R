## modify dRegulation
dRegulation <- function (manifoldOutput) 
{
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl("^X_", geneList)]
  geneList <- gsub("^X_", "", geneList)
  nGenes <- length(geneList)
  eGenes <- nrow(manifoldOutput)/2
  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl("^Y_", eGeneList)]
  eGeneList <- gsub("^Y_", "", eGeneList)
  if (nGenes != eGenes) {
    stop("Number of identified and expected genes are not the same")
  }
  if (!all(eGeneList == geneList)) {
    stop("Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order")
  }
  dMetric <- sapply(seq_len(nGenes), function(G) {
    X <- manifoldOutput[G, ]
    Y <- manifoldOutput[(G + nGenes), ]
    I <- rbind(X, Y)
    O <- dist(I)
    O <- as.numeric(O)
    if (O == 0){
      O <- 1e-17
    }
    return(O)
  })
  lambdaValues <- seq(-2, 2, length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- MASS::boxcox(dMetric ~ 1, plot = FALSE, lambda = lambdaValues)
  BC <- BC$x[which.max(BC$y)]
  if (BC < 0) {
    nD <- 1/(dMetric^BC)
  }
  else {
    nD <- dMetric^BC
  }
  Z <- scale(nD)
  E <- mean(dMetric^2)
  FC <- dMetric^2/E
  pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = "fdr")
  dOut <- data.frame(gene = geneList, distance = dMetric, Z = Z, 
                     FC = FC, p.value = pValues, p.adj = pAdjusted)
  dOut <- dOut[order(dOut$p.value), ]
  dOut <- as.data.frame.array(dOut)
  return(dOut)
}

## fix P-values
# X is manifold alignment results
# d is the dimension of MA needed
# alpha = 0 means both 2d components, alpha = 1 means first, alpha = 2 means back d components
fixPValues <- function(X, d = 2, alpha = 1, gKO = "Trem2"){
  if(alpha == 1){
    X <- X[, 1:d]
  } else if(alpha == 2){
    X <- X[, (ncol(X)/2 + 1):(ncol(X)/2 + d)]
  } else if(alpha == 0){
    X <- X[, c(1:d, (ncol(X)/2 + 1):(ncol(X)/2 + d))]
  }
  X <- as.matrix(X)
  X <- list("manifoldAlignment" = X)
  X$diffRegulation <- dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  if (X$diffRegulation$gene[1] == gKO){
    D <- X$diffRegulation$distance[-1]
    X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
    X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
    X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
    return(X)
  } else{
    return(X)
  }
}

## erichment analysis
enrichFunction <- function(gList, fdrThreshold = 0.05, nCategories = 20){
  library(enrichR)
  E <- enrichr(gList, c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
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
  E <- E[E$Term %in% c('Oxidative phosphorylation','Alzheimer disease','Cholesterol metabolism','Lysosome','Neutrophil mediated immunity'),]
  print(length(unique(E$Term)))
  return(E)
}

## check results
# X is manifold alignment results
# d is the dimension of MA needed
# alpha = 0 means both 2d components, alpha = 1 means first, alpha = 2 means back d components
# threshold for fdr
check_fun <- function(X, d = 2, alpha = 1, fdrThreshold = 0.05, nCategories = 20){
  X <- fixPValues(X, d = d, alpha = alpha)
  # dGenes <- X$diffRegulation$gene[X$diffRegulation$p.value < 0.05]
  dGenes <- X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]
  if (length(dGenes) == 0){
    stop("No significiant genes.")
  }
  # print(c(length(dGenes), length(intersect(markerGenes, dGenes))))
  ## enrichment analysis
  E <- enrichFunction(dGenes, fdrThreshold = fdrThreshold, nCategories = nCategories)
  out_list <- list()
  out_list$gene <- dGenes
  out_list$enrich <- E
  return(out_list)
}

## check sencitivity
check_sensitivity <- function(X, d_index, alpha = 1, fdrThreshold = 0.05, nCategories = 20){
  out_all_list <- list()
  for (i in 1:length(d_index)){
    out_all_list[[i]] <- check_fun(X, d = d_index[i], alpha = alpha, fdrThreshold = fdrThreshold, nCategories = nCategories)
  }
  return(out_all_list)
}

# X is result from check_sencitivity function
plot_sencitivity <- function(X){
  n_list <- length(X)
  n_gene_find <- rep(0, n_list)
  n_path_mean <- rep(0, n_list)
  for (i in 1:n_list){
    tmp <- X[[i]]
    gList <- tmp$gene
    n_gene_find[i] <- length(gList)
    E <- tmp$enrich
    E <- E[E$Term %in% c('Oxidative phosphorylation','Alzheimer disease','Cholesterol metabolism','Lysosome','Neutrophil mediated immunity'),]
    n_path_mean[i] <- length(unique(E$Term))
  }
  A <- ggplot(data = NULL, aes(x = c(1:10), y = n_gene_find)) + 
    geom_line() + geom_point() +
    labs(x = "d", y = "number of genes")
  B <- ggplot(data = NULL, aes(x = c(1:10), y = n_path_mean)) + 
    geom_line() + geom_point() +
    labs(x = "d", y = "number of pathways")
  library(gridExtra)
  return(grid.arrange(A, B, nrow = 1))
}

check_intersect <- function(gList1, gList2){
  return(c(length(gList1), length(gList2), length(intersect(gList1, gList2))))
}

