source("simulation data/utlity.R")

SERGIO <- read.csv('simulation data/SERGIO_create_sim_data/simulationOutput.csv', header = FALSE)
rownames(SERGIO) <- paste0('G', seq_len(nrow(SERGIO)))
colnames(SERGIO) <- paste0('C', seq_len(ncol(SERGIO)))

countMatrix <- SERGIO
set.seed(1)
X <- scTenifoldNet::makeNetworks(countMatrix, q = 0.9)
set.seed(1)
X <- scTenifoldNet::tensorDecomposition(X)
X <- X$X
X <- as.matrix(X)

## not-normalized method
g_res <- fix_pvalue(X = X, gKO = 20, d = 5, alpha = 2, beta = 1) 
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
check_intersect(g_list, g_true)

g_res <- fix_pvalue(X = X, gKO = 50, d = 5, alpha = 2, beta = 1) 
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
check_intersect(g_list, g_true)

g_res <- fix_pvalue(X = X, gKO = 100, d = 5, alpha = 2, beta = 1) 
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
check_intersect(g_list, g_true)

## normalized method
g_res <- fix_pvalue(X = X, gKO = 20, d = 2, alpha = 0, beta = 0, normalize = TRUE) 
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
check_intersect(g_list, g_true)

g_res <- fix_pvalue(X = X, gKO = 41, d = 5, alpha = 0, beta = 0, normalize = TRUE) 
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
check_intersect(g_list, g_true)

g_res <- fix_pvalue(X = X, gKO = 100, d = 2, alpha = 0, beta = 0, normalize = TRUE) 
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
check_intersect(g_list, g_true)

## not-normalized symmetry method
g_res <- fix_pvalue_sym(X = X, gKO = 20, d = 5, alpha = 2) # 1 failed; 0 success.
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
check_intersect(g_list, g_true)

g_res <- fix_pvalue_sym(X = X, gKO = 50, d = 5, alpha = 2) # 1 failed; 0 success.
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
check_intersect(g_list, g_true)

g_res <- fix_pvalue_sym(X = X, gKO = 100, d = 5, alpha = 2) # 1 failed; 0 success.
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
check_intersect(g_list, g_true)

## normalized symmetry method
g_res <- fix_pvalue_sym(X = X, gKO = 20, d = 2, alpha = 2, normalize = TRUE) # 1 failed; 0 success.
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
check_intersect(g_list, g_true)

g_res <- fix_pvalue_sym(X = X, gKO = 50, d = 2, alpha = 2, normalize = TRUE) # 1 failed; 0 success.
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
check_intersect(g_list, g_true)

g_res <- fix_pvalue_sym(X = X, gKO = 100, d = 2, alpha = 2, normalize = TRUE) # 1 failed; 0 success.
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
check_intersect(g_list, g_true)



