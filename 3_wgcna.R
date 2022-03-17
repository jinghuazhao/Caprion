#!/usr/bin/bash

options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(igraph))
suppressMessages(library(RCy3))

load("~/Caprion/pilot/work/es.rda")

wgcna_etc <- function()
# Weighted Correlation Network Analysis
{
  suppressMessages(require(WGCNA))
  enableWGCNAThreads()
# Adjacency matrix using soft thresholding with beta=6
  prot <- t(exprs(protein_all))
  names(prot) <- gsub("_HUMAN","",colnames(prot))
  ADJ <- abs(cor(prot, method="pearson", use="pairwise.complete.obs"))^6
# histogram of k and a scale free topology plot
  k <- as.vector(apply(ADJ,2,sum,na.rm=TRUE))
  sizeGrWindow(10,5)
  par(mfrow=c(1,2))
  hist(k)
  scaleFreePlot(k, main="Check scale free topology\n")
# dissimilarity Topological Overlap Matrix
  dissADJ <- 1 - ADJ
  dissTOM <- TOMdist(ADJ)
  collectGarbage()
# partition around medoids (PAM) based on dissimilarity
  require(cluster)
  for(j in 4:6)
  {
    pam_name <- paste0("pam",j)
    pamTOM_name <- paste0("pamTOM",j)
    assign(pam_name, pam(as.dist(dissADJ),j))
    assign(pamTOM_name,pam(as.dist(dissTOM),j))
    tc <- table(get(pam_name)$clustering,get(pamTOM_name)$clustering)
    print(tc)
    print(diag(tc))
  }
# average linkage hierachical clusterin
# ADJ
  hierADJ <- hclust(as.dist(dissADJ),method="average")
  colorStaticADJ <- as.character(cutreeStaticColor(hierADJ,cutHeight=.99,minSize=5))
  colorDynamicADJ <- labels2colors(cutreeDynamic(hierADJ,method="tree",minClusterSize=5))
  colorDynamicHybridADJ <- labels2colors(cutreeDynamic(hierADJ,distM=dissADJ,cutHeight=0.998,
                                         deepSplit=2,pamRespectsDendro=FALSE))
  colorADJ <- data.frame(pam5$clustering,colorStaticADJ,colorDynamicADJ,colorDynamicHybridADJ)
  pdf("~/Caprion/pilot/work/pamADJ.pdf")
  sizeGrWindow(10,5)
  plotDendroAndColors(dendro=hierADJ,colors=colorADJ,
                      dendroLabels=FALSE,
                      marAll=c(0.2,8,2.7,0.2),
                      main="Gene dendrogram and module colors")
  dev.off()
# TOM
  hierTOM <- hclust(as.dist(dissTOM),method="average");
  colorStaticTOM <- as.character(cutreeStaticColor(hierTOM,cutHeight=.99,minSize=5))
  colorDynamicTOM <- labels2colors(cutreeDynamic(hierTOM,method="tree",minClusterSize=5))
  colorTOM <- data.frame(pamTOM5$clustering,colorStaticTOM,colorDynamicTOM)
  pdf("~/Caprion/pilot/work/pamTOM.pdf")
  plotDendroAndColors(hierTOM,colors=colorTOM,
                      dendroLabels=FALSE,
                      marAll=c(1,8,3,1),
                      main="Gene dendrogram and module colors, TOM dissimilarity")
  dev.off()
  options(width=200)
  cytoscapePing()
  cytoscapeVersionInfo()
  deleteAllNetworks()
  colorADJTOM <- cbind(colorADJ,colorTOM)
  table(colorADJTOM$pamTOM5.clustering)
  for(x in 1:5) print(subset(colorADJTOM,pamTOM5.clustering==x))
  table(colorADJTOM$colorDynamicTOM)
  for(col in c("blue","brown","grey","turquoise","yellow")) print(subset(colorADJTOM,colorDynamicTOM==col))
# Further correlations
  corRaw <- cor(prot,use='pairwise.complete.obs')
  diag(corRaw) <- 0
  distance <- as.dist(1-abs(corRaw))
  colnames(corRaw) <- rownames(corRaw) <- colnames(prot)
  suppressMessages(require(reshape))
  r <- melt(corRaw) %>% mutate(value=ifelse(X1!=X2 & value>=0.7,value,NA))
  colorADJTOM_nogrey <- subset(colorADJTOM,colorStaticTOM!="grey")
  r_nogrey <- melt(corRaw[rownames(colorADJTOM_nogrey),rownames(colorADJTOM_nogrey)]) %>%
              mutate(value=ifelse(X1!=X2 & value>=0.7,value,NA))
  nodes <- data.frame(id=gsub("X4","4",rownames(colorADJTOM_nogrey)),
           group=with(colorADJTOM_nogrey,colorStaticTOM),
           stringsAsFactors=FALSE)
  edges <- data.frame(source=with(r_nogrey,gsub("X4","4",X1)),
           target=with(r_nogrey,gsub("X4","4",X2)),
           weight=with(r_nogrey,value),
           stringsAsFactors=FALSE) %>% filter(!is.na(weight))
  suid_wgnca <- createNetworkFromDataFrames(nodes,edges,title="turquoise", collection="DataFrame")
  layoutNetwork("attribute-circle")
  exportImage("~/Caprion/pilot/work/turquoise.pdf",type="PDF",overwriteFile=TRUE)
  exportNetwork("~/Caprion/pilot/work/turquoise.cyjs","cyjs")
  exportNetwork("~/Caprion/pilot/work/turquoise.sif","SIF")
  exportVisualStyles("~/Caprion/pilot/work/turquoise.json","JSON")
  saveSession("~/Caprion/pilot/work/turquoise.cys")
}

wgcna_etc()
