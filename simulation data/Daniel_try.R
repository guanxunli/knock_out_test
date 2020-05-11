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

#### original method
p_value_KO <- function(gKO){
  Y <- X
  Y[gKO,] <- 0
  MA <- scTenifoldNet::manifoldAlignment(X,Y)
  MA <- MA[, 1:2]
  DR <- dRegulation(MA)
  DR$FC <- DR$distance^2/mean(DR$distance[-seq_len(length(gKO))]^2)
  DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
  DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
  DR <- DR[DR$p.value < 0.05, ]
  return(DR[order(DR$p.value, decreasing = FALSE), ])
}

## Daniel method
g_res <- p_value_KO(20)
g_list <- g_res$gene
g_true <- paste0("G", 16:40)
check_intersect(g_list, g_true)

g_res <- p_value_KO(50)
g_list <- g_res$gene
g_true <- paste0("G", 41:80)
check_intersect(g_list, g_true)

g_res <- p_value_KO(100)
g_list <- g_res$gene
g_true <- paste0("G", 81:100)
check_intersect(g_list, g_true)