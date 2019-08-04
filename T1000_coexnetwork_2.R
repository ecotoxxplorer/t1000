##################################################
## ToxGPrio for Gene Proiritization in Toxicogenomics
## Description: Gene selection, ranking, visualization
## Author: Othman Soufan, othman.soufan@mcgill.ca
## McGill University, Canada. 18 January 2018
## License: GNU GPL (>= 2)
###################################################

require("coexnet")

setwd("~/git/gene-prioritization/geneprioritization/R/T1000/")

TGPMapGenesHSA2Rat <- function(humangenes, ratgenes){
  # This function is to handle mapping especially when Human gene maps to several ones.
  #  It will then make sure to choose the Rat gene that exist in our set.
  #
  # Args:
  #   humangenes: selected genes from Human to map to Rat.
  #   ratgenes: ALL set of genes in a Rat gene expression mat.
  # Returns:
  #   mappedgenes: mapped set of genes
  #
  biomart.gene.mappings <- read_delim("Rat/BioMart_gene_mappings.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
  mappedgenes <- vector(length = length(humangenes))
  hsamapgenes <- vector(length = length(humangenes))
  for(i in c(1:length(humangenes))){
    g <- humangenes[i]
    idx <- which(unlist(biomart.gene.mappings[,2]) == g)
    if(length(idx) == 1){
      if(biomart.gene.mappings[idx,4] %in% ratgenes){
        mappedgenes[i] <- biomart.gene.mappings[idx, 4]
        hsamapgenes[i] <- g
      }
    } else if(length(idx) > 1){
      for(x in idx){
        if(biomart.gene.mappings[x,4] %in% ratgenes){
          mappedgenes[i] <- biomart.gene.mappings[x, 4]
          hsamapgenes[i] <- g
          break;
        }
      }
    }
  }
  mappedgenes <- unlist(mappedgenes)
  hsamapgenes <- unlist(hsamapgenes)
  hsamapgenes[mappedgenes == FALSE] <- NA
  mappedgenes[mappedgenes == FALSE] <- NA
  
  res <- list(mappedgenes, hsamapgenes)
  
  return(res);
}

TGReadData <- function(exprfname, dosefname, org="hsa"){
  require(readr)

  expr.mat <- read_delim(exprfname, " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  expr.genes <- unlist(expr.mat[,1])
    
  if(org == "hsa"){
    humangenes <<- expr.genes
    expr.mat <- data.matrix(expr.mat[,c(2:ncol(expr.mat))])
    rownames(expr.mat) <- expr.genes
    dose.level <- suppressMessages(unlist(read_csv(dosefname, col_names = FALSE)))
    colnames(expr.mat) <- dose.level
  } else if(org == "rat"){
    # Mapping orthologues
    res <- TGPMapGenesHSA2Rat(humangenes, expr.genes)
    mappedgenes <- res[[1]]
    hsamappedgenes <- res[[2]]
    hsamappedgenes <- hsamappedgenes[!is.na(mappedgenes)]
    mappedgenes <- mappedgenes[!is.na(mappedgenes)]
    idx <- match(mappedgenes, expr.genes)
    expr.mat <- expr.mat[idx, ]
    expr.mat[expr.mat == -10000] <- NA
    for(i in 1:nrow(expr.mat)){
      expr.mat[i, is.na(expr.mat[i,])] <- mean(unlist(expr.mat[i,c(2:ncol(expr.mat))]), na.rm = TRUE)
    }
    expr.genes <- hsamappedgenes;#unlist(expr.mat[,1])
    expr.mat <- data.matrix(expr.mat[,c(2:ncol(expr.mat))])
    rownames(expr.mat) <- expr.genes
    dose.level <- suppressMessages(unlist(read_csv(dosefname, col_names = FALSE)))
    colnames(expr.mat) <- dose.level
  }
  
  t <-rep(0,length(dose.level))
  t[dose.level!= "Control"] <- 1
  
  res <- list(Data = expr.mat, t = t)
  
  return(res)
}



# Building the human network
st.time <- Sys.time()
exprfname <- "Human/fullmat_imputed_mean_missingvals"
dosefname <- "Human/dose_level"
res.hsa <- TGReadData(exprfname, dosefname)
set.seed(1)
#sam <-difExprs(expData = res.hsa$Data, treatment = res.hsa$t, fdr = 0.05, DifferentialMethod = "sam")
sam <- res.hsa$Data
cor_pearson_human <- abs(cor(t(sam)))

# # Building the rat (in vitro) network
exprfname <- "Rat/in_vitro/fullmat"
dosefname <- "Rat/in_vitro/dose_level"
humangenes <<- unlist(rownames(sam))
res.rat <- TGReadData(exprfname, dosefname, "rat")
cor_pearson_rat_invitro <- abs(cor(t(res.rat$Data)))
# 
# Building the rat (in vivo) network
exprfname <- "Rat/in_vivo/Liver/Single/fullmat"
dosefname <- "Rat/in_vivo/Liver/Single/dose_level"
res.rat.vivo <- TGReadData(exprfname, dosefname, "rat")
cor_pearson_rat_invivo <- abs(cor(t(res.rat.vivo$Data)))

common.genes <- intersect(intersect(rownames(cor_pearson_human), rownames(cor_pearson_rat_invitro)), rownames(cor_pearson_rat_invivo))

cor_pearson_human <- cor_pearson_human[match(common.genes, rownames(cor_pearson_human)),match(common.genes, rownames(cor_pearson_human))]
cor_pearson_rat_invitro <- cor_pearson_rat_invitro[match(common.genes, rownames(cor_pearson_rat_invitro)),match(common.genes, rownames(cor_pearson_rat_invitro))]
cor_pearson_rat_invivo <- cor_pearson_rat_invivo[match(common.genes, rownames(cor_pearson_rat_invivo)),match(common.genes, rownames(cor_pearson_rat_invivo))]

cor.thresh <- 0.5
diag(cor_pearson_human) <- 0
cor_pearson_human.tmp <- cor_pearson_human
cor_pearson_human[cor_pearson_human < cor.thresh] <- 0
cor_pearson_human[cor_pearson_human >= cor.thresh] <- 1

diag(cor_pearson_rat_invitro) <- 0
cor_pearson_rat_invitro.tmp <- cor_pearson_rat_invitro
cor_pearson_rat_invitro[cor_pearson_rat_invitro < cor.thresh] <- 0
cor_pearson_rat_invitro[cor_pearson_rat_invitro >= cor.thresh] <- 1

diag(cor_pearson_rat_invivo) <- 0
cor_pearson_rat_invivo.tmp <- cor_pearson_rat_invivo
cor_pearson_rat_invivo[cor_pearson_rat_invivo < cor.thresh] <- 0
cor_pearson_rat_invivo[cor_pearson_rat_invivo >= cor.thresh] <- 1

cor.final <- cor_pearson_human+cor_pearson_rat_invitro+cor_pearson_rat_invivo
cor.final.tmp <- cor_pearson_human.tmp+cor_pearson_rat_invitro.tmp+cor_pearson_rat_invivo.tmp
cor.final[is.na(cor.final)] <- 0
cor.final.tmp[is.na(cor.final.tmp)] <- 0

cor.final.tmp[cor.final != 3] <- 0

cor.final.tmp.v1 <- colSums(cor.final.tmp)
cor.final.tmp <- cor.final.tmp[cor.final.tmp.v1!=0, cor.final.tmp.v1!=0]
genes.rank <- common.genes[cor.final.tmp.v1!=0]
cor.final.tmp.v1 <- cor.final.tmp.v1[cor.final.tmp.v1!=0]
o <- order(cor.final.tmp.v1, decreasing = TRUE)

genes.rank <- common.genes[o]


###################################
# Prepare prior scores
###################################
ranked.genes.ctd <- unlist(read_csv("~/git/gene-prioritization/geneprioritization/R/T1000/ranked_genes_ctd.txt", 
                             col_names = FALSE))
ranked.genes.ctd.priors <- c(length(ranked.genes.ctd):1)
ranked.genes.ctd.priors <- ranked.genes.ctd.priors/max(ranked.genes.ctd.priors)
all.priors <- numeric(length = nrow(cor.final.tmp))+0.00001

priors.idx <- match(rownames(cor.final.tmp), ranked.genes.ctd)
priors.idx <- priors.idx[!is.na(priors.idx)]
cor.final.idx <- match(ranked.genes.ctd, rownames(cor.final.tmp))
cor.final.idx <- cor.final.idx[!is.na(cor.final.idx)]

all.priors[cor.final.idx] <- ranked.genes.ctd.priors[priors.idx]

###################################
# Construct igraph object for clustering
###################################
cor.final.tmp2 <- cor.final.tmp
#cor.final.tmp2[lower.tri(cor.final.tmp2)] <- 0
cor.final.tmp2 <- cor.final.tmp2/max(cor.final.tmp2)

idx.p <- which(lower.tri(cor.final.tmp2), arr.ind = TRUE)

for(i in 1:nrow(idx.p)){
  gene_a <- idx.p[i,1]
  gene_b <- idx.p[i,2]
  prior <- (all.priors[gene_a]+all.priors[gene_b])/2
  cor.final.tmp2[gene_a, gene_b] <- cor.final.tmp2[gene_a, gene_b]*prior
  cor.final.tmp2[gene_b, gene_a] <- cor.final.tmp2[gene_b, gene_a]*prior
}

# # Prepare files for running MCL on the server
# idx.p <- which(cor.final.tmp2 != 0, arr.ind = TRUE)
# mcl.abs <- matrix(, nrow = nrow(idx.p), ncol=3)
# gg <- rownames(cor.final.tmp2)
# for(i in 1:nrow(idx.p)){
#   gene_a <- idx.p[i,1]
#   gene_b <- idx.p[i,2]
#   mcl.abs[i,] <- c(gg[gene_a], gg[gene_b], cor.final.tmp2[gene_a, gene_b])
# }
# mcl.abs <- as.data.frame(mcl.abs)
# write.table(mcl.abs, "GraphClustering/MCL/bin/coex_tggates.abc", row.names = FALSE, col.names = FALSE, quote = FALSE)

require('MCL')
st.mcl <- Sys.time(); res.cluster <- mcl(x = cor.final.tmp2, addLoops = TRUE); print(Sys.time()-st.mcl)

# require(igraph)
# rowcols <- which(cor.final.tmp2>0, arr.ind = TRUE)
# g <- cbind(rownames(cor.final.tmp2)[rowcols[,1]], rownames(cor.final.tmp2)[rowcols[,2]])
# weights <- cor.final.tmp2[rowcols]
# 
# graph <- graph_from_data_frame(g, directed = FALSE)
# E(graph)$weight <- weights
# 
# comms <- cluster_louvain(graph)

Sys.time()-st.time