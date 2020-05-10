## fix P-values
# X is manifold alignment results
# d is the dimension of MA needed
# alpha = 0 means both 2d components, alpha = 1 means first, alpha = 2 means back d components
fixPValues <- function(X, d = 2, alpha = 1){
  if(alpha == 1){
    X <- X[, 1:d]
  } else if(alpha == 2){
    X <- X[, (ncol(X)/2 + 1):(ncol(X)/2 + d)]
  } else if(alpha == 0){
    X <- X[, c(1:d, (ncol(X)/2 + 1):(ncol(X)/2 + d))]
  }
  X <- list("manifoldAlignment" = X)
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-1]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
  return(X)
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
  ## select Set
  # selectedSet <- rep(FALSE, nrow(E))
  # for(i in seq_len(nrow(E))){
  #   if(i == 1){
  #     selectedSet[i] <- TRUE
  #   } else {
  #     A <- unique(unlist(strsplit(E[which(selectedSet[seq_len(i)]),'Genes'], ';')))
  #     B <- unlist(strsplit(E[i,'Genes'], ';'))
  #     selectedSet[i] <- !all(B %in% A)
  #   }
  # }
  # gSets <- table(toupper(E$Term))
  # gSets <- names(gSets[gSets > 1])
  # for(i in gSets){
  #   selectedSet[which(toupper(E$Term) %in% i)[-1]] <- FALSE
  # }
  # E <- E[selectedSet,]
  # if(nrow(E) > nCategories){
  #   E <- E[seq_len(nCategories),]  
  # }
  
  ## print results
  name_need <- c("GO_Cellular_Component_2018.1", "Reactome_2016.1", "BioPlanet_2019.1",
                 "GO_Cellular_Component_2018.4", "GO_Cellular_Component_2018.2")
  term_need <- c("Lamellar body", "Surfactant metabolism", "Cell adhesion molecules", "Multivesicular body lumen", "Lysosome")
  index1 <- which(rownames(E) %in% name_need)
  index2 <- which(E$Term %in% term_need)
  index <- intersect(index1, index2)
  print(c(nrow(E), index))
  return(E)
}

## check results
# X is manifold alignment results
# d is the dimension of MA needed
# alpha = 0 means both 2d components, alpha = 1 means first, alpha = 2 means back d components
# threshold for fdr
check_fun <- function(X, d = 2, alpha = 1, fdrThreshold = 0.05, nCategories = 20){
  X <- fixPValues(X, d = d, alpha = alpha)
  dGenes <- X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]
  # print(c(length(dGenes), length(intersect(markerGenes, dGenes))))
  ## enrichment analysis
  E <- enrichFunction(dGenes, fdrThreshold = fdrThreshold, nCategories = nCategories)
  out_list <- list()
  out_list$gene <- dGenes
  out_list$enrich <- E
  return(out_list)
}

check_sensitivity <- function(X, d_index, alpha = 1, fdrThreshold = 0.05, nCategories = 20){
  out_all_list <- list()
  for (i in 1:length(d_index)){
    out_all_list[[i]] <- check_fun(X, d = d_index[i], alpha = alpha, fdrThreshold = fdrThreshold, nCategories = nCategories)
  }
  return(out_all_list)
}

check_intersect <- function(gList1, gList2){
  return(c(length(gList1), length(gList2), length(intersect(gList1, gList2))))
}