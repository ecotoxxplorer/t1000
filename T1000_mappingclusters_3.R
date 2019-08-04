##################################################
## ToxGPrio for Gene Proiritization in Toxicogenomics
## Description: Gene selection, ranking, visualization
## Author: Othman Soufan, othman.soufan@mcgill.ca
## McGill University, Canada. 18 January 2018
## License: GNU GPL (>= 2)
###################################################

# map output from MCL algorithm
library(readr)
clusters <- read_delim("~/git/gene-prioritization/geneprioritization/R/T1000/GraphClustering/out.coex_tggates_2.abc.I33", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)
View(clusters)

cluster.sel <- c(6:8)
aftercluster2 <- aftercluster
for(i in c(1:length(cluster.sel))){
  cidx <- cluster.sel[i]
  cgenes <- unlist(clusters[cidx,])
  cgenes <- cgenes[!is.na(cgenes)]
  for(j in c((i+1):length(cluster.sel))){
    xidx <- cluster.sel[j]
    xgenes <- unlist(clusters[xidx,])
    xgenes <- xgenes[!is.na(xgenes)]
    
    for(g1 in cgenes){
      for(g2 in xgenes){
        idx1 <- which(rownames(aftercluster2) == g1)
        idx2 <- which(rownames(aftercluster2) == g2)
        aftercluster2[idx1, idx2] <- 0
      }
    }
  }
}
