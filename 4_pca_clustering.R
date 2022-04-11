options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(dplyr))
suppressMessages(library(gap))

load("~/Caprion/pilot/work/es.rda")
prot <- t(exprs(protein_all))
colnames(prot) <- gsub("_HUMAN","",colnames(prot))

ZYQ.na <- c("BROX","CT027","GHRL","PSB6")
UDP.na <- c("BROX","NCF2","SEM7A")

digitalnames <- colnames(prot)[grepl("^[0-9]",colnames(prot))]
cat("digital names:",digitalnames,"\n")
cat("X-prefixed names:",gsub("(^[0-9])","X\\1",digitalnames),"\n")

pca_clustering <- function(prot)
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
  list(pca=pca,pc1pc2=pc1pc2,eigenvec=eigenvec,km=km,mc=mc)
}

pca_clustering_plot <- function(pca,mc,out)
{
  suppressMessages(library(plotly))
  trace <- list(
    mode = "markers", 
    name = "f", 
    type = "scatter3d", x = with(pca,x[,2]),y = with(pca,x[,1]), z = with(pca,x[,3]),
    frame = NULL, 
    marker = list(line = list(color = "transparent"), color = "rgba(102,194,165,1)", fillcolor = "rgba(102,194,165,0.5)")
  )
  layout <- list(
    scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")), 
    xaxis = list(domain = c(0, 1)), 
    yaxis = list(domain = c(0, 1)), 
    margin = list(b = 40, l = 60, r = 10, t = 25),
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

  htmlwidgets::saveWidget(p,file=out)
}

quantro_sparsenetgls <- function(edata,batch,mod)
{
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
  if (FALSE) {
     y <- t(edata)
     sngls <- sparsenetgls(y,mod)
     save(sngls,file="~/Caprion/pilot/work/sparsenetgls.rda")
     plotsngls(sngls,ith_lambda=5)
  }
}

normalise_lr <- function(d,batches)
{
  normfun <- function(col,verbose=FALSE)
  {
    if (verbose) cat(names(d[col]),col,"\n")
    y <- invnormal(d[[col]])
    l <- lm(y~mod,data=d[mnames])
    r <- y-predict(l,na.action=na.pass)
    invnormal(r)
  }
  if (batches==1) {
     mod <- model.matrix(as.formula(paste0(c("~agePulse","sexPulse",paste0("PC",1:20)),collapse="+")), data=d)
     covars <- c(names(d)[grep("agePulse|sexPulse",names(d))],paste0("PC",1:20))
     proteins <- setdiff(names(d),c("FID","IID","batch",paste0("ppc",1:3),covars))
  } else {
     mod <- model.matrix(as.formula(paste0(c("~agePulse","sexPulse",paste0("ppc",1:3),paste0("PC",1:20)),collapse="+")), data=d)
     covars <- c(names(d)[grep("agePulse|sexPulse",names(d))],paste0("ppc",1:3),paste0("PC",1:20))
     proteins <- setdiff(names(d),c("FID","IID","batch",covars))
  }
  mnames <- colnames(mod)[-1]
  z <- sapply(names(d[proteins]), normfun)
  colnames(z) <- names(d[proteins])
  rownames(z) <- d[["IID"]]
  caprion_lr <- data.frame(d[c("FID","IID")],z)
  names(caprion_lr) <- gsub("^X([0-9])","\\1",names(caprion_lr))
  write.table(caprion_lr,file=paste0("~/Caprion/pilot/work/caprion-",batches,".pheno"),quote=FALSE,row.names=FALSE,sep="\t")
}

normalise <- function(prot)
{
  pca_km_mc <- pca_clustering(prot)
  attach(pca_km_mc)
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
     dev.off()
     ZYQ_mc <- read.csv("~/Caprion/pilot/ZYQ_PC1_groups_20200703.csv")
     mc_ZYQ_mc <- cbind(ZYQ_mc,classification=with(mc,classification)[grepl("ZYQ",names(mc$classification))])
     with(mc_ZYQ_mc,table(pc1_group,classification))
     pca_clustering_plot(pca,mc,"~/Caprion/pilot/work/pca_clustering.html")
  })
# Phenotype files
  data <- read.csv("~/Caprion/INTERVALdata_15SEP2021.csv")
  eigenvec <- read.delim("~/Caprion/pilot/data/merged_imputation.eigenvec")
  pilotsMap <- read.csv("~/Caprion/pilotsMap_15SEP2021.csv")
  OmicsMap <- read.csv("~/Caprion/INTERVAL_OmicsMap_20210915.csv")
  grouping <- data.frame(caprion_id=names(with(mc,classification)),classification=with(mc,classification)) %>%
              left_join(data.frame(with(pca,x)[,1:3],caprion_id=rownames(with(pca,x)))) %>%
              rename(ppc1=PC1,ppc2=PC2,ppc3=PC3)
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
  detach(pca_km_mc)
  write.table(select(pheno,Affymetrix_gwasQC_bl),file="~/Caprion/pilot/work/caprion.id",quote=FALSE,row.names=FALSE,col.names=FALSE)
  caprion_pheno <- mutate(pheno,FID=Affymetrix_gwasQC_bl,IID=Affymetrix_gwasQC_bl) %>%
                   select(FID,IID,gsub("(^[0-9])","X\\1",colnames(prot)))
  names(caprion_pheno) <- gsub("^X([0-9])","\\1",names(caprion_pheno))
  write.table(caprion_pheno,file="~/Caprion/pilot/work/caprion.pheno",quote=FALSE,row.names=FALSE,sep="\t")

  suppressMessages(library(sva))
  mod <- model.matrix(as.formula(paste0(c("~agePulse","sexPulse","ppc1","ppc2","ppc3",paste0("PC",1:20)),collapse="+")), data=pheno)
  mnames <- colnames(mod)[-1]
  mdat <- pheno[mnames]
  mcol <- apply(mdat,2,is.na)
  many <- !apply(mcol,1,any)

# 1. parametric adjustment
  combat_edata1 <- ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)

# 2. non-parametric adjustment, mean-only version
  combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

  edata <- edata[,many]
  batch <- batch[many]
  quantro_sparsenetgls(edata,batch,mod)
# 3. reference-batch version, with covariates
  combat_edata3 <- ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3, prior.plots=TRUE)

# invnormal transformation
  allvars <- c("FID","IID","agePulse","sexPulse","batch",paste0("ppc",1:3),paste0("PC",1:20),gsub("(^[0-9])","X\\1",colnames(prot)))
  dat <- mutate(pheno,FID=Affymetrix_gwasQC_bl,IID=Affymetrix_gwasQC_bl)[allvars]
  for(batches in 1:3)
  {
    d <- filter(dat[many,],batch==batches)
    normalise_lr(d,batches)
    write.table(d[c("FID","IID")],file=paste0("~/Caprion/pilot/work/caprion-",batches,".id"),col.names=FALSE,row.names=FALSE)
  }
# return normalised data from above
  list(edata=t(combat_edata3),batch=batch)
}

edata_batch <- normalise(prot)
attach(edata_batch)
combat_quantro <- quantro(t(edata),batch,B=10000)
quantroPlot(combat_quantro)
pca_km_mc <- pca_clustering(edata)
detach(edata_batch)
attach(pca_km_mc)
with(mc, if (interactive()) {rgl::plot3d(with(pca,x[,c(2,1,3)]),col=classification)})
pca_clustering_plot(pca,mc,"~/Caprion/pilot/work/pca_clustering_combat.html")
detach(pca_km_mc)
