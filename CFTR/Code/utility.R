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
  X <- as.matrix(X)
  X <- list("manifoldAlignment" = X)
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment)
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  D <- X$diffRegulation$distance[-1]
  X$diffRegulation$FC <- X$diffRegulation$distance^2/mean(D^2)
  X$diffRegulation$p.value <- pchisq(X$diffRegulation$FC, df = 1, lower.tail = FALSE)
  X$diffRegulation$p.adj <- p.adjust(X$diffRegulation$p.value, method = 'fdr')
  return(X)
}

## check results
# X is manifold alignment results
# d is the dimension of MA needed
# alpha = 0 means both 2d components, alpha = 1 means first, alpha = 2 means back d components
# threshold for fdr
check_fun <- function(X, d = 2, alpha = 1, fdrThreshold = 0.05, nCategories = 20){
  X <- fixPValues(X, d = d, alpha = alpha)
  dGenes <- X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]
  out_list <- list()
  out_list$gene <- dGenes
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
  for (i in 1:n_list){
    tmp <- X[[i]]
    gList <- tmp$gene
    n_gene_find[i] <- length(gList)
  }
  A <- ggplot(data = NULL, aes(x = c(1:10), y = n_gene_find)) + 
    geom_line() + geom_point() +
    labs(x = "d", y = "number of genes")
  return(A)
}

check_intersect <- function(gList1, gList2){
  return(c(length(gList1), length(gList2), length(intersect(gList1, gList2))))
}

