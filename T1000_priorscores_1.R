##################################################
## ToxGPrio for Gene Proiritization in Toxicogenomics
## Description: Gene selection, ranking, visualization
## Author: Othman Soufan, othman.soufan@mcgill.ca
## McGill University, Canada. 18 January 2018
## License: GNU GPL (>= 2)
###################################################

# Params: 
#   symbolslist: global variable of list of genes (e.g. c("ZFP57", "TGM3", "SLC39A2", "AMFR", "NOL4", "DOCK9", "MON2"))
TGPPerformEnrichAnalysis <- function(type,deg.vec = symbolslist){
  require(qusage)
  
  if(type == "kegg"){
    file.gmt <- "KEGG_gene_sets.gmt"
  }
  else if(type == "reactome"){
    file.gmt <- "Reactome_gene_sets.gmt"
  }
  else if(type == "gobp"){
    file.gmt <- "GO_biological_process.gmt"
  }
  else if(type == "gomf"){
    file.gmt <- "GO_molecular_function.gmt"
  }
  else if(type == "gocc"){
    file.gmt <- "GO_cellular_component.gmt"
  }
  else if(type == "cancer"){
    file.gmt <- "cancer_gene_neighborhoods.gmt"
  }
  else if(type == "chemgenetic"){
    file.gmt <- "chemical_and_genetic_perturbations.gmt"
  }
  else if(type == "immuno"){
    file.gmt <- "immunologic_signatures.gmt"
  }
  else if(type == "hallmark"){
    file.gmt <- "hallmark_gene_sets.gmt"
  }
  
  anot.type <- strsplit(file.gmt, ".gmt")[[1]][1]
  file.gmt <- paste0("MSigDB/", file.gmt)
  anot.set <- read.gmt(file.gmt)
  
  # prepare for the result table
  set.size<-length(anot.set);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(anot.set);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes
  anot.universe <- unique(unlist(anot.set))
  hits.inx <- deg.vec %in% anot.universe;
  deg.vec <- deg.vec[hits.inx];
  
  q.size<-length(deg.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(anot.set, 
                       function(x) {
                         deg.vec[deg.vec%in%unlist(x)];
                       }
  );
  
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
  
  # total unique gene number
  uniq.count <- length(anot.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(anot.set, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] <= 0.05;
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }
  
  #get gene symbols
  enrichres <- as.data.frame(res.mat)
  enrichres$Terms <- rownames(enrichres)
  res <- c(enrichres$Terms[1], enrichres$P.Value[1])
  return(res)
  
  #file.nm <- paste0(anot.type,".txt")
  #write.table(data.frame(Terms=rownames(res.mat), res.mat),file=file.nm,sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  #return(1);
}


TGPInitDataObjects <- function(){
  # Constructs a global dataSet and genes objects for storing data .
  #
  # Args:
  #   None.
  #
  # Returns:
  #   The sample covariance between x and y.
  #
  # TODO: Consider data initialization for other species (e.g. rat)
  
  require(igraph)
  require(ggfortify)
  require(factoextra)
  
  if(!(exists("dataSet") && exists("genes"))){
    # Load CTD pre-built graph with human and rat subgraphs extracted
    # Therapeutic relationships are excluded from original graph and only
    # mechanistic interactions that possibily indicate toxicities are kept.
    load("Data/toxchem_gene_graph.RData")
    # Load KEGG pathways and Hallmark genesets for mapping to input list of genes
    load("Data/pathways_data.RData")
    
    # Prepare graph of genes-chems-disease(toxic endpoints)
    edges.chem.gene <- as.vector(t(data.frame(human.chemgene.graph$X2, human.chemgene.graph$X4)))
    edges.chem.disease <- as.vector(t(data.frame(ctd_chem_graph$X1, ctd_chem_graph$X2)))
    edges <- c(edges.chem.gene, edges.chem.disease)
    graph <- make_graph(edges, directed = FALSE)
    
    # Work with human genes
    genes <<- unique(human.chemgene.graph$X4)
    
    # Feature extraction - phase I: define genes context in the toxicogenomics space
    # Compute closeness of genes in the CTD graph
    vals.closeness <- closeness(graph)
    vals.names <- names(vals.closeness)
    genes.closeness <- unlist(lapply(genes, function(x) vals.closeness[which(vals.names == x)]))
    # Compute betweeness of genes in the CTD graph
    vals.betweeness <- betweenness(graph, directed = FALSE)
    vals.names <- names(vals.betweeness)
    genes.betweeness <- unlist(lapply(genes, function(x) vals.betweeness[which(vals.names == x)]))
    # Compute degree of genes in the CTD graph
    vals.degree <- degree(graph)
    vals.names <- names(vals.degree)
    genes.degree <- unlist(lapply(genes, function(x) vals.degree[which(vals.names == x)]))
    
    # Feature extraction - phase II: define genes context in the biological space
    # Prepare pathways features
    genes.pathways <- matrix(0,length(genes), nrow(hallmark.genesets)+nrow(kegg.pathways))
    for(i in c(1:length(genes))){
      idx <- which(hallmark.genesets == genes[i], arr.ind = TRUE)
      genes.pathways[i,idx[,1]] <- 1
      idx <- which(kegg.pathways == genes[i], arr.ind = TRUE)
      genes.pathways[i,(idx[,1]+nrow(hallmark.genesets))] <- 1
    }
    
    genes.betweeness[genes.betweeness == 0] <- 1
    # Scale and assign data
    dataSet.orig <<- data.frame(genes.closeness, genes.betweeness, genes.degree, genes.pathways)
    dataSet <<- data.frame(log(genes.closeness), log(genes.betweeness), log(genes.degree), scale(genes.pathways))
    
  }
}

TGPDrawFigure <- function(allres, k=3, outsz=10){
  tiff(filename = "Results/Figures/ToxGPrio_filtering_steps.tiff",
       width = 178, height = 170, units = "mm", res = 300)
  genes.sel <- integer()
  data <- dataSet
  allidx <- c(1:nrow(data))
  layout(matrix(c(1,1,1,2,2,2,3,4,5), 3,3, byrow = TRUE))#par(mfrow=c(3,1))
  i <- 0
  set.seed(1)
  
  iteration <- 2
  
  while(i < iteration){
    # PCA Analysis
    res.pca <- prcomp(data)
    res.ind <- get_pca_ind(res.pca)
    
    # Clustering Analysis
    res.pca.data <- data.frame(res.ind$coord[,1], res.ind$coord[,2])
    kmeans.result <- kmeans(res.pca.data, centers=k)
    centers <- kmeans.result$centers[kmeans.result$cluster, ] # "centers" is a data frame of 3 centers but the length of iris dataset so we can canlculate distance difference easily.
    distances <- sqrt(rowSums((res.pca.data - centers)^2))
    outliers <- order(distances, decreasing=T)[1:outsz]
    centroids <- order(distances)[1:outsz]
    
    # Meta-data for outliers
    outliers.CTD <- colMeans(dataSet.orig[allidx[outliers], 1:3])[3]
    outliers.hallmark <- mean(colMeans(dataSet.orig[allidx[outliers], 4:53]))*100
    outliers.kegg <- mean(colMeans(dataSet.orig[allidx[outliers], 54:239]))*100
    # Meta-data for centroids
    centroids.CTD <- colMeans(dataSet.orig[allidx[centroids], 1:3])[3]
    centroids.hallmark <- mean(colMeans(dataSet.orig[allidx[centroids], 4:53]))*100
    centroids.kegg <- mean(colMeans(dataSet.orig[allidx[centroids], 54:239]))*100
    
    # For highly correlated genes, keep one of them
    corscores <- cor(t(data[outliers,]))
    corscores[lower.tri(corscores, diag = TRUE)] <- 0
    simidx <- which(corscores > 0.89, arr.ind = TRUE)
    if(length(simidx) > 0){
      genes.sel <- c(genes.sel, allidx[outliers[-c(simidx[,1])]])
    } else{
      genes.sel <- c(genes.sel, allidx[outliers])
    }
    
    # Retrieve outliers as interesting set of genes
    genes.sel <- c(genes.sel, outliers)
    
    data <- data[-outliers,]
    allidx <- allidx[-outliers]
    
    title <- paste0("Itr", i+1)
    colors <- kmeans.result$cluster
    platte1 <- rgb(253/255, 117/255, 111/255)
    platte2 <- rgb(0/255, 186/255, 69/255)
    platte3 <- rgb(93/255, 153/255, 251/255)
    colors[colors == 1] <- platte1
    colors[colors == 2] <- platte2
    colors[colors == 3] <- platte3
    
    plot(res.pca.data, pch=19, col=colors, cex=1, xlab="PC1", ylab="PC2", main="", cex.lab=1.5, cex.axis=1.5)
    title(main = list(title, cex = 1,
                      col = "black", font = 2))
    #colors <- c(platte1, platte2, platte3)
    points(kmeans.result$centers, col="black", pch=15, cex=1.6)
    points(res.pca.data[outliers,], pch="+", col=4, cex=3)
    if(i == 0){
      text(res.pca.data[outliers[3],1]+5, res.pca.data[outliers[3],2]+50, paste0("Gene degree (CTD): ", outliers.CTD), font=2)
      text(res.pca.data[outliers[3],1]+5, res.pca.data[outliers[3],2]+38, paste0("Hallmark: ", round(outliers.hallmark, 2),"%"), font=2)
      text(res.pca.data[outliers[3],1]+5, res.pca.data[outliers[3],2]+23, paste0("KEGG: ", round(outliers.kegg, 2), "%"), font=2)
      
      text(res.pca.data[centroids[2],1]-10, res.pca.data[centroids[2],2]-30, paste0("Gene degree (CTD): ", centroids.CTD), font=2)
      text(res.pca.data[centroids[2],1]-10, res.pca.data[centroids[2],2]-43, paste0("Hallmark: ", round(centroids.hallmark, 2),"%"), font=2)
      text(res.pca.data[centroids[2],1]-10, res.pca.data[centroids[2],2]-58, paste0("KEGG: ", round(centroids.kegg, 2), "%"), font=2)
    } else if(i == 1){
      text(res.pca.data[outliers[3],1]-3, res.pca.data[outliers[3],2]+53, paste0("Gene degree (CTD): ", outliers.CTD), font=2)
      text(res.pca.data[outliers[3],1]-3, res.pca.data[outliers[3],2]+41, paste0("Hallmark: ", round(outliers.hallmark, 2),"%"), font=2)
      text(res.pca.data[outliers[3],1]-3, res.pca.data[outliers[3],2]+26, paste0("KEGG: ", round(outliers.kegg, 2), "%"), font=2)
      
      text(res.pca.data[centroids[2],1]+3, res.pca.data[centroids[2],2]-10, paste0("Gene degree (CTD): ", centroids.CTD), font=2)
      text(res.pca.data[centroids[2],1]+3, res.pca.data[centroids[2],2]-23, paste0("Hallmark: ", round(centroids.hallmark, 2),"%"), font=2)
      text(res.pca.data[centroids[2],1]+3, res.pca.data[centroids[2],2]-38, paste0("KEGG: ", round(centroids.kegg, 2), "%"), font=2)
    }
    i <- i + 1
  }
  plot(log(allres[[2]][,1]), pch=19, col=platte3, cex=1, xlab="Iteration", ylab="Gene degree (CTD)", main="", cex.lab=1.5, cex.axis=1.5, type="o")
  lines(log(allres[[2]][,4]),col=platte1)
  legend(50,max(log(allres[[2]][,1])), c("Out.","Cen."), lwd=c(1.5,1.5),col=c(platte3, platte1))
  plot(allres[[2]][,2], pch=19, col=platte3, cex=1, xlab="Iteration", ylab="Hallmark Matching %", main="", cex.lab=1.5, cex.axis=1.5, type="o")
  lines(allres[[2]][,5],col=platte1)
  legend(50,max(allres[[2]][,2]), c("Out.","Cen."), lwd=c(1.5,1.5),col=c(platte3, platte1))
  plot(allres[[2]][,3], pch=19, col=platte3, cex=1, xlab="Iteration", ylab="KEGG Pathway Matching %", main="", cex.lab=1.5, cex.axis=1.5, type="o", lty=c(1,1), lwd=c(2.5,2.5))
  lines(allres[[2]][,6],col=platte1)
  legend(50,max(allres[[2]][,3]), c("Out.","Cen."), lwd=c(1.5,1.5),col=c(platte3, platte1))
  
  dev.off()
  
  return(genes.sel)
}

TGPFilteringGenes <- function(expr.mat, dose.level, k=3, outsz=10, thresh=0.25){
  genes.sel <- integer()
  data <- dataSet
  allidx <- c(1:nrow(data))
  set.seed(1)
  corsz <- 0
  
  # Initialization of gene scoring variables
  genes.all <- rownames(expr.mat)
  genes.scores <- list()
  
  datasz <- nrow(data) # Original data size
  while(nrow(data) > 140){
    print(paste0("Number of genes ranked ", length(genes.sel)))
    # PCA Analysis
    res.pca <- prcomp(data)
    res.ind <- get_pca_ind(res.pca)
    
    # Clustering Analysis
    res.pca.data <- data.frame(res.ind$coord[,1], res.ind$coord[,2])
    kmeans.result <- kmeans(res.pca.data, centers=k)
    centers <- kmeans.result$centers[kmeans.result$cluster, ] # "centers" is a data frame of 3 centers but the length of iris dataset so we can canlculate distance difference easily.
    distances <- sqrt(rowSums((res.pca.data - centers)^2))
    outliers <- order(distances, decreasing=T)[1:outsz]
    centroids <- order(distances)[1:outsz]

    # For highly correlated genes, keep one of them
    # corscores <- cor(t(data[outliers,]))
    # corscores[lower.tri(corscores, diag = TRUE)] <- 0
    # simidx <- which(corscores > 0.89, arr.ind = TRUE)
    # if(length(simidx) > 0){
    #   corsz <- corsz + length(unique(c(simidx[,1])))
    #   idxset <- allidx[outliers[-c(simidx[,1])]]
    # } else{
    #   idxset <- allidx[outliers]
    # }
    idxset <- allidx[outliers]# Comment if to use previous code for correlation
    genes.sel <- c(genes.sel, idxset)
    
    # Map cluster information to gene expression space and run Limma
    genenm.idx <- match(genes[idxset], rownames(expr.mat))
    genenm.idx <- genenm.idx[!is.na(genenm.idx)]
    if(sum(!is.na(genenm.idx) == TRUE) > 0){
      sub.expr <- expr.mat[genenm.idx, ]
      if(is.vector(sub.expr)){
        topgenek <- 1
        limma.res <- TGPLimmaFeaSel(sub.expr, dose.level, "Treatment", "Control", topgenek)
        sub.idx <- which(genes.all == row.names(expr.mat)[genenm.idx])
        genes.scores[genes.all[sub.idx]] <- limma.res$P.Value
      } else{
        topgenek <- nrow(sub.expr)
        limma.res <- TGPLimmaFeaSel(sub.expr, dose.level, "Treatment", "Control", topgenek)
        sub.idx <- as.vector(sapply(row.names(limma.res), function(x) which(genes.all == x)))
        genes.scores[genes.all[sub.idx]] <- limma.res$P.Value
      }
    }
    
    data <- data[-outliers,]
    allidx <- allidx[-outliers]
  }
  print(corsz)
  allres <- list(Genes = genes.sel, Scores = genes.scores)
  
  return(allres)
}

TGPLimmaFeaSel <- function(exprmat, dose.level, trt.dose.level, ctr.dose.level, topgenek){
  require(limma)
  
  # Condition preparation
  cond <- rep(0,length(dose.level));
  
  # initialize
  if(trt.dose.level == "Treatment"){
    # initialize
    t1 <- cond; c <- cond; 
    t1[dose.level!="Control"] <- 1;
    c[dose.level==ctr.dose.level] <- 1;
    design=cbind(treatment=t1, control=c)
    cont.matrix <- makeContrasts(STvsCO=treatment-control, levels=design)
  }
  else if(trt.dose.level == "High vs. Low"){
    # initialize
    t1 <- cond; t2 <- cond; c <- cond;
    t1[dose.level=="High"] <- 1;
    if("Low" %in% dose.level){# Some files do not have tests for Low doses
      t2[dose.level=="Low"] <- 1;
    }
    else if("Middle" %in% dose.level){
      t2[dose.level=="Middle"] <- 1;
    }
    else{
      return -1;
    }
    c[dose.level==ctr.dose.level] <- 1;
    design=cbind(high=t1, low=t2, control=c)
    cont.matrix <- makeContrasts(high-low,high-control,low-control, levels=design)
  }
  else{ # For only 2 groups of any type of treatment compared to another single type of control
    # initialize
    t1 <- cond; 
    t1[dose.level==trt.dose.level] <- 1;
    c <- cond; c[dose.level==ctr.dose.level] <- 1;
    design=cbind(treatment=t1, control=c)
    cont.matrix <- makeContrasts(STvsCO=treatment-control, levels=design)
  }
  
  # limma fitting
  fit <- lmFit(exprmat, design);
  fit <- contrasts.fit(fit, cont.matrix);
  efit <- eBayes(fit);
  
  # Results
  res <- topTable(efit,number=topgenek, adjust.method = "fdr");
  #selfeaidx <- as.integer(rownames(res))
  
  return(res);
}

TGPPrepMultiOmics <- function(expr.mat, genes.human, max.ldh=104, min.ldh=96){
  # Parse and impute missing values of Rat in vivo gene expression matrix
  rat.invivo.expr <- read_delim("Rat/in_vivo/Liver/Single/fullmat", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  genes.rat <- TGPMapGenesHSA2Rat(genes.human, unlist(rat.invivo.expr[,1]))
  idx <- match(genes.rat, unlist(rat.invivo.expr[,1]))
  rat.invivo.expr <- rat.invivo.expr[idx,]
  rat.invivo.expr[rat.invivo.expr == -10000] <- NA
  for(i in 1:nrow(rat.invivo.expr)){
    rat.invivo.expr[i, is.na(rat.invivo.expr[i,])] <- mean(unlist(rat.invivo.expr[i,c(2:ncol(rat.invivo.expr))]), na.rm = TRUE)
  }
  
  # Parse and impute missing values of Rat in vitro gene expression matrix
  rat.invitro.expr <- read_delim("Rat/in_vitro/fullmat", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  genes.rat <- TGPMapGenesHSA2Rat(genes.human, unlist(rat.invitro.expr[,1]))
  idx <- match(genes.rat, unlist(rat.invitro.expr[,1]))
  rat.invitro.expr <- rat.invitro.expr[idx,]
  rat.invitro.expr[rat.invitro.expr == -10000] <- NA
  for(i in 1:nrow(rat.invitro.expr)){
    rat.invitro.expr[i, is.na(rat.invitro.expr[i,])] <- mean(unlist(rat.invitro.expr[i,c(2:ncol(rat.invitro.expr))]), na.rm = TRUE)
  }
  
  require(omicade4)
  # Get LDH survival rate and does reads for all experiments (i.e. Human, Rat (in vitro and in vivo))
  #  Annotate levels of LDH as toxic based on input low and high ldh rates
  ldh.human <- unlist(read_csv("Human/y", col_names = FALSE))
  doses.human <- unlist(read_csv("Human/dose", col_names = FALSE))
  #genes.human <<- unlist(expr.mat[,1])
  expr.mat <- expr.mat[,c(2:ncol(expr.mat))]
  # Exclude experiments with missing/NA reported LDH levels
  expr.mat <- expr.mat[,!is.na(ldh.human)]
  ldh.human <- ldh.human[!is.na(ldh.human)]
  labels.human <- ldh.human
  labels.human[ldh.human > max.ldh | ldh.human < min.ldh] <- "Toxic"
  labels.human[labels.human!="Toxic"] <- "Non-Toxic"
  
  # Rat (in vitro)
  ldh.rats <- unlist(read_csv("Rat/in_vitro/y", col_names = FALSE))
  doses.rats <- unlist(read_csv("Rat/in_vitro/dose", col_names = FALSE))
  genes.rat.invitro <<- unlist(rat.invitro.expr[,1])
  genes.rat.invitro <- as.vector(unlist(biomart.gene.mappings[match(genes.rat.invitro, unlist(biomart.gene.mappings[,4])),2]))
  rat.invitro.expr <- rat.invitro.expr[,c(2:ncol(rat.invitro.expr))]
  rat.invitro.expr <- rat.invitro.expr[,!is.na(ldh.rats)]
  ldh.rats <- ldh.rats[!is.na(ldh.rats)]
  labels.rat.invitro <- ldh.rats
  labels.rat.invitro[ldh.rats > max.ldh | ldh.rats < min.ldh] <- "Toxic"
  labels.rat.invitro[labels.rat.invitro!="Toxic"] <- "Non-Toxic"
  
  # Rat (in vivo)
  genes.rat.invivo <<- unlist(rat.invivo.expr[,1])
  genes.rat.invivo <- as.vector(unlist(biomart.gene.mappings[match(genes.rat.invivo, unlist(biomart.gene.mappings[,4])),2]))
  rat.invivo.expr <- rat.invivo.expr[,c(2:ncol(rat.invivo.expr))]
  pathology.rats <- unlist(read_csv("Rat/in_vivo/Liver/Single/y", col_names = FALSE))
  doses.rats.invivo <- unlist(read_csv("Rat/in_vivo/Liver/Single/dose", col_names = FALSE))
  labels.rat.invivo <- pathology.rats
  labels.rat.invivo[pathology.rats == 2] <- "Toxic"
  labels.rat.invivo[pathology.rats != 2] <- "Non-Toxic"
  
  sz <- 156
  sz2 <- 300
  idx.toxic <- which(labels.human == "Toxic")
  idx.nontoxic <- which(labels.human == "Non-Toxic")
  human.momics <- data.frame(expr.mat[,idx.toxic[1:sz]], expr.mat[,idx.nontoxic[1:sz2]])
  rownames(human.momics) <- genes.human
  
  idx.toxic <- which(labels.rat.invitro == "Toxic")
  idx.nontoxic <- which(labels.rat.invitro == "Non-Toxic")
  rat.invitro.momics <- data.frame(rat.invitro.expr[,idx.toxic[1:sz]], rat.invitro.expr[,idx.nontoxic[1:sz2]])
  rownames(rat.invitro.momics) <- genes.human
  
  idx.toxic <- which(labels.rat.invivo == "Toxic")
  idx.nontoxic <- which(labels.rat.invivo == "Non-Toxic")
  rat.invivo.momics <- data.frame(rat.invivo.expr[,idx.toxic[1:sz]], rat.invivo.expr[,idx.nontoxic[1:sz2]])
  rownames(rat.invivo.momics) <- genes.human
  
  tggates.momics <- list(Human=human.momics, "Rat in vitro" = rat.invitro.momics, "Rat in vivo" = rat.invivo.momics)
  
  
  phenovec <- rep("Toxic", sz+sz2)
  phenovec[(sz+1):(sz+sz2)] <- "Non-Toxic"
  
  omicsdata <- list(data = tggates.momics, phenovec = phenovec)
  
  return(omicsdata)
}

TGPMomicsPrioritize <- function(omicsdata, mcoin, txgenesz=15){
  sz1 <- nrow(omicsdata$data$Human)
  sz2 <- nrow(omicsdata$data$`Rat in vitro`)
  sz3 <- nrow(omicsdata$data$`Rat in vivo`)
  Tco.a <- mcoin$mcoa$Tco[(1:sz1),];
  Tco.b <- mcoin$mcoa$Tco[(sz1+1):(sz1+sz2),];
  Tco.c <- mcoin$mcoa$Tco[(sz1+sz2+1):(sz1+sz2+sz3),]
  
  # a1.lim=c(-Inf, 0) specifies directionality of toxicity in sample space. At the moment, the expert specifies it.
  #  TODO: set a1.lim in automated fashion based on region where toxic samples cluster
  toxicity.genes <- selectVar(mcoin, a1.lim=c(-Inf, 0), a2.lim=c(-Inf, Inf))
  # TODO: choose txgenesz for top genes contributing to toxicity in an automated fashion.
  #  One can choose top k genes contributing to cluster of genes in the toxicity region of sample space.
  gene.stat <- plotVar(mcoin, var=toxicity.genes$var[1:txgenesz], var.lab=TRUE)
  selgenes <- as.vector(gene.stat$Variables)
  
  genes.a <- rownames(omicsdata$data$Human)
  avg.a <- vector(length = 10)
  for(g in selgenes){
    idx <- which(genes.a == g)
    avg.a <- avg.a + Tco.a[idx,]
  }
  avg.a <- avg.a/length(selgenes)
  scores.a <- vector(length = length(genes.a))
  for(i in c(1:length(genes.a))){
    scores.a[i] <- 1/sqrt(rowSums((Tco.a[i,] - avg.a)^2))
  }
  
  genes.b <- rownames(omicsdata$data$`Rat in vitro`)
  avg.b <- vector(length = 10)
  for(g in selgenes){
    idx <- which(genes.b == g)
    if(length(idx) != 0){
      avg.b <- avg.b + Tco.b[idx,]
    }
    else{
      avg.b <- avg.b + 0
    }
  }
  avg.b <- avg.b/length(selgenes)
  scores.b <- vector(length = length(genes.b))
  for(i in c(1:length(genes.b))){
    scores.b[i] <- 1/sqrt(rowSums((Tco.b[i,] - avg.b)^2))
  }
  
  genes.c <- rownames(omicsdata$data$`Rat in vivo`)
  avg.c <- vector(length = 10)
  for(g in selgenes){
    idx <- which(genes.c == g)
    if(length(idx) != 0){
      avg.c <- avg.c + Tco.b[idx,]
    }
    else{
      avg.c <- avg.c + 0
    }
  }
  avg.c <- avg.c/length(selgenes)
  scores.c <- vector(length = length(genes.c))
  for(i in c(1:length(genes.c))){
    scores.c[i] <- 1/sqrt(rowSums((Tco.c[i,] - avg.c)^2))
  }
  
  # Assume rat genes are always part of human
  #  Also, assume gene set for rat are the same for in vitro and in vivo
  for(g in genes.a){
    if(g %in% genes.b){
      idx1 <- which(genes.a == g)
      idx2 <- which(genes.b == g)
      idx3 <- which(genes.c == g)
      scores.a[idx1] <- (scores.a[idx1]+scores.b[idx2]+scores.c[idx3])/3
    } else{
      idx1 <- which(genes.a == g)
      scores.a[idx1] <- (scores.a[idx1])
    }
  }
  
  scores <- scores.a#(scores.a+scores.b+scores.c)/3
  
  arrange.idx <- order(scores, decreasing=T)
  arrange.scores <- scores[arrange.idx]
  arrange.genes <- genes.a[arrange.idx]
  
  return(arrange.genes)
}

convertHumanToRat <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x ,
                   mart = human, attributesL = c("rgd_symbol"), martL = rat, uniqueRows=T)
  
  output <- unique(genesV2[, 2])
  
  return(output)
}

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
  biomart.gene.mappings <<- read_delim("Rat/BioMart_gene_mappings.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
  mappedgenes <- vector(length = length(humangenes))
  for(i in c(1:length(humangenes))){
    g <- humangenes[i]
    idx <- which(unlist(biomart.gene.mappings[,2]) == g)
    if(length(idx) == 1){
      if(biomart.gene.mappings[idx,4] %in% ratgenes){
        mappedgenes[i] <- biomart.gene.mappings[idx, 4]
      }
    } else if(length(idx) > 1){
      for(x in idx){
        if(biomart.gene.mappings[x,4] %in% ratgenes){
          mappedgenes[i] <- biomart.gene.mappings[x, 4]
          break;
        }
      }
    }
  }
  mappedgenes <- unlist(mappedgenes)
  mappedgenes[mappedgenes == FALSE] <- NA
  return(mappedgenes);
}

#setwd("~/git/gene-prioritization/geneprioritization/R/ToxGPrio/")
#TGPInitDataObjects()
#require(readr)
#exprfname <- "Human/fullmat_imputed_mean_missingvals"
#expr.mat <- read_delim(exprfname, " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
#expr.genes <- unlist(expr.mat[,1])
#expr.mat <- data.matrix(expr.mat[,c(2:ncol(expr.mat))])
#rownames(expr.mat) <- expr.genes
#dosefname <- "Human/dose_level"
#dose.level <- suppressMessages(unlist(read_csv(dosefname, col_names = FALSE)))
#colnames(expr.mat) <- dose.level

# The following 5 lines are the ones to start running !!!!!!!!!!!!!!!!!
# allres <- TGPFilteringGenes(expr.mat, dose.level)
# ranked.genes <- genes[allres$Genes]
# write.table(ranked.genes, "ranked_genes_ctd.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# pvals <- unlist(allres[[2]])
# write.table(pvals, "ranked_genes_ctd_pvals.txt", row.names = TRUE, col.names = FALSE, quote = FALSE)

#sc <- allres[[2]]
#genes.all <- rownames(expr.mat)
#b <- order(sc)
#topgenes <- genes.all[b]
#TGPDrawFigure(allres)

# # Step 3: Next mapp genes to gene expression dataset and apply Limma
# genes.sel <- genes[allres[[1]]]

# idx <- match(genes.sel, unlist(expr.mat[,1]))
# genes.sel <- genes.sel[!is.na(idx)]
# idx <- idx[!is.na(idx)]
# expr.mat <- expr.mat[idx, ]
# dosefname <- "Human/dose_level"
# dose.level <- suppressMessages(read_csv(dosefname, col_names = FALSE))
# topgenek <- 384
# trt.dose.level <- "Treatment"
# ctr.dose.level <- "Control"
# genes.limma <- TGPLimmaFeaSel(expr.mat[,c(2:ncol(expr.mat))], dose.level, trt.dose.level, ctr.dose.level, topgenek)
# expr.mat <- expr.mat[genes.limma,]
# genes.human <- genes.sel[genes.limma]

# Step 4: Ranking based on multi-omics over final selected 384 genes for ranking
# Prepare multi-omics datasets
#omicsdata <- TGPPrepMultiOmics(expr.mat, genes.human)
#mcoin <- mcia(omicsdata$data, cia.nf=10)
#plot(mcoin, axes=1:2, phenovec=omicsdata$phenovec, sample.lab=FALSE, df.color=1:3)

# Prioritize genes based on contribution of all experiments towards toxicity
#genes.ranked <- TGPMomicsPrioritize(omicsdata, mcoin)

