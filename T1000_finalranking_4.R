# setwd("~/git/gene-prioritization/geneprioritization/R/T1000/")
# 
# 
# library(readr)
# genes.names <- unlist(read_csv("~/git/gene-prioritization/geneprioritization/R/T1000/genes_names_2.txt", 
#                         col_names = FALSE))
# genes.scores <- unlist(read_csv("~/git/gene-prioritization/geneprioritization/R/T1000/genes_scores_2.txt", 
#                                col_names = FALSE))
# 
# 
# mcl.clusters <- as.matrix(read_delim("~/git/gene-prioritization/geneprioritization/R/T1000/GraphClustering/out.coex_tggates_2.abc.I33", 
#                            "\t", escape_double = FALSE, col_names = FALSE, 
#                            trim_ws = TRUE))
# 
# 
# priors.pvals <- read_delim("~/git/gene-prioritization/geneprioritization/R/T1000/ranked_genes_ctd_pvals.txt", 
#                                      " ", escape_double = FALSE, col_names = FALSE, 
#                                      trim_ws = TRUE)
# 
# 
# genes.priors.orig <- as.matrix(read_delim("~/git/gene-prioritization/geneprioritization/R/T1000/genes_priors_original_data.txt", 
#                                          " ", escape_double = FALSE, col_names = FALSE, 
#                                          trim_ws = TRUE))
# 
# 
# genes.cor.sets <- read_delim("~/git/gene-prioritization/geneprioritization/R/T1000/correlated_gene_sets_prior.txt", 
#                                          " ", escape_double = FALSE, col_names = FALSE, 
#                                          trim_ws = TRUE)
# 
# mcl.clusters.tmp <- mcl.clusters
# for(i in c(1:nrow(mcl.clusters.tmp))){
#   clus <- mcl.clusters.tmp[i,]
#   for(j in c(1:nrow(genes.cor.sets))){
#     target <- unlist(genes.cor.sets[j, ])
#     target <- target[!is.na(target)]
#     idx <- match(target, clus)
#     idx <- idx[!is.na(idx)]
#     if(length(idx) > 1){
#       clus[idx[2:length(idx)]] <- NA
#     }
#   }
#   mcl.clusters.tmp[i, ] <- clus
# }
# 
# write.table(mcl.clusters.tmp, "out.coex_tggates_2.abc.I33_updated", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
# 
# ###################Final feature ranking
library(readr)
genes_rankv3 <- read_csv("~/git/gene-prioritization/geneprioritization/R/T1000/Results3/genes_names_4.txt",
                         col_names = FALSE)
genes_rankv3 <- unlist(genes_rankv3)

genes_scoresv3 <- read_csv("~/git/gene-prioritization/geneprioritization/R/T1000/Results3/genes_scores_4.txt",
                           col_names = FALSE)
genes_scoresv3 <- unlist(genes_scoresv3)

cluster_perfv3 <- read_csv("~/git/gene-prioritization/geneprioritization/R/T1000/Results3/cluster_perf_4.txt",
                           col_names = FALSE)
cluster_perfv3 <- unlist(cluster_perfv3)


c.sz <- rowSums(!is.na(mcl.clusters2))
c.sz.tmp <- c.sz
idx <- order(cluster_perfv3, decreasing = TRUE)
genes.final <- c()
tmp <- genes.final
g.overall.sz <- 1500
overall.sz <- g.overall.sz
genes_scoresv3.tmp <- genes_scoresv3
while(length(genes.final) < g.overall.sz){
  cluster.feasz <- ceiling((cluster_perfv3/sum(cluster_perfv3))*overall.sz)
  
  for(i in c(1:length(idx))){
    s <- idx[i]
    if(s == 1){
      st <- 1
      end <- c.sz[(idx[i])]
    } else{
      st <- sum(c.sz[1:(idx[i]-1)])+1
      end <- st+c.sz[(idx[i])]-1
    }
    
    print(paste0(i, " ", st, " ", end))
    
    l.genes.rank <- order(genes_scoresv3.tmp[st:end], decreasing = TRUE)
    all.idx <- st:end
    all.idx <- all.idx[l.genes.rank]
    l.genes.names <- genes_rankv3[st:end]
    l.genes.names <- l.genes.names[l.genes.rank]
    
    if(sum(genes_scoresv3.tmp[st:end]) != 0){
      if(cluster.feasz[idx[i]] > 0 && c.sz.tmp[idx[i]] > 0){
        if(cluster.feasz[idx[i]] > c.sz.tmp[idx[i]]){
          genes.final <- c(genes.final, as.vector(l.genes.names[1:c.sz.tmp[idx[i]]]))
          genes_scoresv3.tmp[all.idx[1:c.sz.tmp[idx[i]]]] <- 0
          c.sz.tmp[idx[i]] <- c.sz.tmp[idx[i]]-c.sz.tmp[idx[i]]
        } else {
          genes.final <- c(genes.final, as.vector(l.genes.names[1:cluster.feasz[idx[i]]]))
          genes_scoresv3.tmp[all.idx[1:cluster.feasz[idx[i]]]] <- 0
          c.sz.tmp[idx[i]] <- c.sz.tmp[idx[i]]-cluster.feasz[idx[i]]
        }
      }
    }
    
  }
  
  overall.sz <- g.overall.sz-length(genes.final)
  print(length(genes.final))
}

# # Step1: For each cluster, remove duplicate genes based on the prior
# mcl.clusters.updated <- list()
# for(c in mcl.clusters){
#   selidx <- match(c, genes.priors.orig[,1])
#   selidx <- selidx[!is.na(selidx)]
#   genes <- genes.priors.orig[selidx, 1] 
#   data <- genes.priors.orig[selidx, 2:ncol(genes.priors.orig)]
#   data <- data.matrix(as.data.frame(data))
#   data.cor <- cor(t(data))
#   data.cor[lower.tri(data.cor, diag = TRUE)] <- 0
#   simidx <- which(data.cor > 0.99, arr.ind = TRUE)
#   if(length(simidx) > 0){
#     genes.sel <- c(genes.sel, allidx[outliers[-c(simidx[,1])]])
#   } else{
#     genes.sel <- c(genes.sel, allidx[outliers])
#   }
#   
# }