options(width=200)
suppressMessages(library(Biobase))

load("~/Caprion/pilot/work/es.rda")
prot <- t(exprs(protein_all))
colnames(prot) <- gsub("_HUMAN","",colnames(prot))

ZYQ.na <- c("BROX","CT027","GHRL","PSB6")
UDP.na <- c("BROX","NCF2","SEM7A")

pca_clustering_plot <- function(pca,mc)
{
  suppressMessages(library(plotly))
  trace <- list(
    mode = "markers", 
    name = "f", 
    type = "scatter3d", x = with(pca,x[,2]),y = with(pca,x[,1]), z = with(pca,x[,3]),
    frame = NULL, 
    marker = list(
      line = list(color = "transparent"), 
      color = "rgba(102,194,165,1)", 
      fillcolor = "rgba(102,194,165,0.5)"
    )
  )
  layout <- list(
    scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")), 
    xaxis = list(domain = c(0, 1)), 
    yaxis = list(domain = c(0, 1)), 
    margin = list(
      b = 40, 
      l = 60, 
      r = 10, 
      t = 25
    ), 
    hovermode = "closest", 
    showlegend = TRUE
  )

  suppressMessages(library(ggplot2))
  p <- plot_ly() %>%
       add_trace(mode=trace$mode, name=trace$name, type=trace$type, x=trace$x, y=trace$y, z=trace$z, frame=trace$frame, marker=trace$marker) %>%
       layout(scene=layout$scene, xaxis=layout$xaxis, yaxis=layout$yaxis, margin=layout$margin, hovermode=layout$hovermode, showlegend=layout$showlegend)

  p <- with(list(x=with(pca,x[,2]),y=with(pca,x[,1]),z=with(pca,x[,3]),color=mc$classification,id=with(pca,rownames(x))),
       plot_ly(x = ~x, y = ~y, z = ~z, color = ~color, colors=c("black","red")) %>%
       add_markers(type="scatter3d", text=id) %>%
       layout(scene = list(xaxis = list(title = "PC2", tickmode = "array", autotick = TRUE, tickfont = list(size = 10)), 
                           yaxis = list(title = "PC1", tickmode = "array", autotick = TRUE, tickfont = list(size = 10)),
                           zaxis = list(title = "PC3", tickfont = list(size = 10)),
                           aspectratio = list(x = 0.9, y = 1, z = 0.6)),
                           title = "3d plot",
                           showlegend = TRUE))

  htmlwidgets::saveWidget(p,file="~/Caprion/pilot/work/pca_clustering.html")
}

pca_clustering <- function()
{
  all <- prot[,!colnames(prot)%in%union(ZYQ.na,UDP.na)]
  pca <- prcomp(all, rank=50, scale=TRUE)
  pc1pc2 <- with(pca,x)[,1:2]
  rownames(pc1pc2) <- rownames(all)
  eigenvec <- with(pca,rotation)[,1:2]
  cor(eigenvec)
  cor(pc1pc2)
  pdf("~/Caprion/pilot/work/pca_clustering.pdf")
  screeplot(pca, npcs=20, type="lines", main="PCA screeplot")
  plot(eigenvec,pch=19,cex=0.6)
  title("Eigenvectors")
  plot(pc1pc2,pch=19,cex=0.6)
  title("Principal components")
  biplot(pca,cex=0.1)
  title("biplot")
# K-means clustering
  km <- kmeans(pc1pc2,2)
  table(with(km,cluster))
  plot(pc1pc2, col = with(km,cluster), pch=19, cex=0.8)
  points(with(km,centers), col = 1:2, pch = 8, cex = 2)
  title("K-means clustering")
# Model-based clustering
  library(mclust)
  mc <- Mclust(pc1pc2,G=2)
  summary(mc)
  table(with(mc,classification))
  plot(mc, what=c("classification"))
  title("Model-based clustering")
  ZYQ_mc <- read.csv("ZYQ_PC1_groups_20200703.csv")
  mc_ZYQ_mc <- cbind(ZYQ_mc,classification=with(mc,classification)[grepl("ZYQ",names(mc$classification))])
  with(mc_ZYQ_mc,table(pc1_group,classification))
  if (interactive()) with(mc,
  {
     png(file.path("~/Caprion/pilot/work/3d.png"), res=300, width=12, height=10, units="in")
     scatterplot3d::scatterplot3d(with(pca,x[,c(2,1,3)]), color=c("blue","red")[classification], main="Plot of the PC1, PC2 and PC3", pch=16)
     legend("right", legend=levels(as.factor(classification)), col=c("blue", "red"), pch=16)
     dev.off()
     rgl::plot3d(with(pca,x[,c(2,1,3)]),col=classification)
  })
  dev.off()
  pca_clustering_plot(pca,mc)
# Phenotype files
  pilotsMap <- read.csv("pilotsMap_15SEP2021.csv")
  OmicsMap <- read.csv("INTERVAL_OmicsMap_20210915.csv")
  data <- read.csv("INTERVALdata_15SEP2021.csv")
  id <- c("identifier","Affymetrix_gwasQC_bl","caprion_id")
  date <- c("attendanceDate","sexPulse","monthPulse","yearPulse","agePulse")
  covars <- c("ethnicPulse","ht_bl","wt_bl","CRP_bl","TRANSF_bl","processDate_bl","processTime_bl","classification")
  grouping <- data.frame(caprion_id=names(with(mc,classification)),classification=with(mc,classification))
  id_date_covars <- merge(merge(data,merge(pilotsMap,OmicsMap,by="identifier",all=TRUE),by="identifier",all=TRUE),grouping,by="caprion_id")
  dim(id_date_covars)
  head(id_date_covars[c(id,date,covars)])
}

pca_clustering()
