options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(dplyr))

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
# K-means clustering
  km <- kmeans(pc1pc2,2)
  table(with(km,cluster))
# Model-based clustering
  library(mclust)
  mc <- Mclust(pc1pc2,G=2)
  summary(mc)
  table(with(mc,classification))
  with(mc, if (interactive()) {rgl::plot3d(with(pca,x[,c(2,1,3)]),col=classification)} else {
     pdf("~/Caprion/pilot/work/pca_clustering.pdf")
     scatterplot3d::scatterplot3d(with(pca,x[,c(1,2,3)]), color=c("blue","red")[classification], main="Plot of the PC1, PC2 and PC3", pch=16)
     legend("right", legend=levels(as.factor(classification)), col=c("blue", "red"), pch=16)
     screeplot(pca, npcs=20, type="lines", main="PCA screeplot")
     plot(eigenvec,pch=19,cex=0.6)
     title("Eigenvectors")
     plot(pc1pc2,pch=19,cex=0.6)
     title("Principal components")
     biplot(pca,cex=0.1)
     title("biplot")
     plot(pc1pc2, col = with(km,cluster), pch=19, cex=0.8)
     points(with(km,centers), col = 1:2, pch = 8, cex = 2)
     title("K-means clustering")
     plot(mc, what=c("classification"))
     title("Model-based clustering")
     ZYQ_mc <- read.csv("~/Caprion/ZYQ_PC1_groups_20200703.csv")
     mc_ZYQ_mc <- cbind(ZYQ_mc,classification=with(mc,classification)[grepl("ZYQ",names(mc$classification))])
     with(mc_ZYQ_mc,table(pc1_group,classification))
     dev.off()
     pca_clustering_plot(pca,mc)
  })
# Phenotype files
  data <- read.csv("~/Caprion/INTERVALdata_15SEP2021.csv")
  eigenvec <- read.delim("~/Caprion/pilot/data/merged_imputation.eigenvec")
  pilotsMap <- read.csv("~/Caprion/pilotsMap_15SEP2021.csv")
  OmicsMap <- read.csv("~/Caprion/INTERVAL_OmicsMap_20210915.csv")
  grouping <- data.frame(caprion_id=names(with(mc,classification)),classification=with(mc,classification)) %>%
              left_join(data.frame(pc1pc2,caprion_id=rownames(pc1pc2))) %>%
              rename(ppc1=PC1,ppc2=PC2)
  id <- c("identifier","Affymetrix_gwasQC_bl","caprion_id")
  date <- c("attendanceDate","sexPulse","monthPulse","yearPulse","agePulse")
  covars <- c("ethnicPulse","ht_bl","wt_bl","CRP_bl","TRANSF_bl","processDate_bl","processTime_bl","classification")
  id_date_covars <- merge(merge(data,merge(pilotsMap,OmicsMap,by="identifier",all=TRUE),by="identifier",all=TRUE),grouping,by="caprion_id")
  dim(id_date_covars)
  head(id_date_covars[c(id,date,covars)])
  pheno <- select(OmicsMap,identifier,Affymetrix_gwasQC_bl,caprion_id) %>%
           merge(eigenvec,by.x="Affymetrix_gwasQC_bl",by.y="X.FID",all.x=TRUE) %>%
           left_join(grouping) %>%
           left_join(data) %>%
           right_join(data.frame(prot,caprion_id=rownames(prot))) %>%
           mutate(batch=match(substr(sampleNames(protein_all),1,3),c("ZWK","ZYQ","UDP")))
  edata <- exprs(protein_all)
  rownames(edata) <- sub("_HUMAN","",rownames(edata))
  edata <- edata[!rownames(edata)%in%union(ZYQ.na,UDP.na),]
  batch <- pheno$batch

  suppressMessages(library(sva))
  mod <- model.matrix(as.formula(paste0(c("~agePulse","sexPulse","classification",paste0("PC",1:20)),collapse="+")), data=pheno)

# 1. parametric adjustment
  combat_edata1 <- ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)

# 2. non-parametric adjustment, mean-only version
  combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

  n <- 2491
  edata <- edata[,1:n]
  batch <- head(batch,n)
# 3. reference-batch version, with covariates
  combat_edata3 <- ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3, prior.plots=TRUE)

  suppressMessages(library(quantro))
  suppressMessages(library(doParallel))
  registerDoParallel(cores=10)
  quantro(edata,batch,B=10000)
  png("~/Caprion/pilot/work/matboxplot.png",width=12,height=10,unit="in",res=300)
  par(cex=0.6,cex.lab=0.6)
  matboxplot(edata,batch)
  dev.off()
  png("~/Caprion/pilot/work/matdensity.png",width=12,height=10,unit="in",res=300)
  matdensity(edata, batch, xlab = " ", ylab = "density", ylim=c(0,2),
             main = "Protein levels", brewer.n = 8, brewer.name = "Dark2")
  legend('top', c("ZWK", "ZYQ", "UDP"), col = c(1, 2, 3), lty = 1, lwd = 3)
  dev.off()
  suppressMessages(library(sparsenetgls))
  y <- t(edata)
  sngls <- sparsenetgls(y,mod)
  plotsngls(sngls,ith_lambda=5)
}

pca_clustering()
