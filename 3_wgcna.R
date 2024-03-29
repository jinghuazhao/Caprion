options(width=200)
suppressMessages(library(Biobase))

load("~/Caprion/pilot/work/es.rda")
prot <- t(exprs(protein_all))
colnames(prot) <- gsub("_HUMAN","",colnames(prot))

mcol <- apply(prot,2,is.na)
many <- !apply(mcol,1,any)
s <- prot[many,]
r <- cor(s,s)
g <- as.matrix(r>=0.7)+0

suppressMessages(library(RColorBrewer))
suppressMessages(library(heatmaply))
heatmaply(r,k_col=5,k_row=5,colors = colorRampPalette(brewer.pal(3, "RdBu"))(256),file="~/Caprion/pilot/work/corr.html")

suppressMessages(library(ComplexHeatmap))
f1 <- Heatmap(r,row_names_gp = gpar(fontsize = 3),column_names_gp = gpar(fontsize = 3))
f2 <- Heatmap(g,row_names_gp = gpar(fontsize = 3),column_names_gp = gpar(fontsize = 3))
f <- f1
f
ggplot2::ggsave("heatmap.pdf")

wgcna_etc <- function()
# Weighted Correlation Network Analysis
{
  suppressMessages(library(cluster))
  suppressMessages(library(WGCNA))
  enableWGCNAThreads()
# Adjacency matrix using soft thresholding with beta=6
  ADJ <- abs(cor(prot, method="pearson", use="pairwise.complete.obs"))^6
# k <- softConnectivity(prot)
# connectivity is the sum of the adjacency to the other nodes.
# sizeGrWindow(10,5)
# histogram of k and a scale free topology plot
  k <- as.vector(apply(ADJ,2,sum,na.rm=TRUE))
  pdf("~/Caprion/pilot/work/k.pdf",height=8,width=12)
  par(mfrow=c(1,2))
  hist(k)
  scaleFreePlot(k, main="Check scale free topology\n")
  dev.off()
# dissimilarity Topological Overlap Matrix
  dissADJ <- 1 - ADJ
  dissTOM <- TOMdist(ADJ)
  collectGarbage()
# partition around medoids (PAM) based on dissimilarity
  for(j in 4:8)
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
  colorDynamicHybridADJ <- labels2colors(cutreeDynamic(hierADJ,distM=dissADJ,cutHeight=0.998,deepSplit=2,pamRespectsDendro=FALSE))
  colorADJ <- data.frame(pam5$clustering,colorStaticADJ,colorDynamicADJ,colorDynamicHybridADJ)
  pdf("~/Caprion/pilot/work/pamADJ.pdf")
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
  colorADJTOM <- cbind(colorADJ,colorTOM)
  table(colorADJTOM$pamTOM6.clustering)
  for(x in 1:5) print(subset(colorADJTOM,pamTOM5.clustering==x))
  table(colorADJTOM$colorDynamicTOM)
  Colors <- c("blue","brown","grey","turquoise","yellow")
  for(col in Colors) print(subset(colorADJTOM,colorDynamicTOM==col))
# A reframe
  pwr <- c(1:10, seq(from=12, to=30, by=2))
  sft <- pickSoftThreshold(prot, powerVector=pwr, verbose=5)
  ADJ <- abs(cor(prot, method="pearson", use="pairwise.complete.obs"))^6
  dissADJ <- 1-ADJ
  dissTOM <- TOMdist(ADJ)
  TOM <- TOMsimilarityFromExpr(prot)
  Tree <- hclust(as.dist(1-TOM), method="average")
  for(j in pwr)
  {
    pam_name <- paste0("pam",j)
    assign(pam_name, pam(as.dist(dissADJ),j))
    pamTOM_name <- paste0("pamTOM",j)
    assign(pamTOM_name,pam(as.dist(dissTOM),j))
    tc <- table(get(pam_name)$clustering,get(pamTOM_name)$clustering)
    print(tc)
    print(diag(tc))
  }
  colorStaticTOM <- as.character(cutreeStaticColor(Tree,cutHeight=.99,minSize=5))
  colorDynamicTOM <- labels2colors(cutreeDynamic(Tree,method="tree",minClusterSize=5))
  Colors <- data.frame(pamTOM5$clustering,colorStaticTOM,colorDynamicTOM)
  pdf("~/Caprion/pilot/work/reframe.pdf")
  plotDendroAndColors(Tree, Colors, dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
  dev.off()
  meg <- moduleEigengenes(prot, color=1:ncol(prot), softPower=6)
}

wgcna_etc()

cytoscape <- function()
# Cytoscape
{
  suppressMessages(library(cluster))
  suppressMessages(library(dplyr))
  suppressMessages(library(RCy3))
  suppressMessages(library(reshape))
  cytoscapePing()
  cytoscapeVersionInfo()
  deleteAllNetworks()
  corRaw <- cor(prot,use='pairwise.complete.obs')
  distance <- as.dist(1-abs(corRaw))
  colnames(corRaw) <- rownames(corRaw) <- colnames(prot)
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
  getLayoutNames()
  layoutNetwork("isom")
  exportImage("turquoise.pdf",type="PDF",overwriteFile=TRUE)
  exportNetwork("turquoise.cyjs","cyjs")
  exportNetwork("turquoise.sif","SIF")
  exportVisualStyles("turquoise.json","JSON")
  saveSession("turquoise.cys")
}

# https://www.datanovia.com/en/blog/how-to-create-a-beautiful-interactive-heatmap-in-r/
