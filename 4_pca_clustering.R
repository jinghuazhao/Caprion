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
  p <- plot_ly()
  p <- with(trace,add_trace(p,mode=mode, name=name, type=type, x=x, y=y, z=z, frame=frame, marker=marker))
  p <- with(layout,layout(p,scene=scene, xaxis=xaxis, yaxis=yaxis, margin=margin, hovermode=hovermode, showlegend=showlegend))
  p2 <- with(list(x=with(pca,x[,2]),y=with(pca,x[,1]),z=with(pca,x[,3]),color=mc$classification,id=with(pca,rownames(x))),
        plot_ly(x = ~x, y = ~y, z = ~z, color = ~color, colors=c("black","red")) %>%
        add_markers(type="scatter3d", text=id) %>%
        layout(scene = list(xaxis = list(title = "PC2", tickmode = "array", autotick = TRUE, tickfont = list(size = 10)), 
                            yaxis = list(title = "PC1", tickmode = "array", autotick = TRUE, tickfont = list(size = 10)),
                            zaxis = list(title = "PC3", tickfont = list(size = 10)),
                            aspectratio = list(x = 0.9, y = 1, z = 0.6)),
                            title = "3d plot",
                            showlegend = TRUE))
  if (interactive()) rgl::plot3d(with(pca,x[,c(2,1,3)]),col=with(mc,classification)) else htmlwidgets::saveWidget(p2,file=out)
}

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
  if(!interactive()) {
     pdf(paste0("~/Caprion/analysis/output/pca_clustering",suffix,".pdf"))
     scatterplot3d::scatterplot3d(with(pca,x[,c(1,2,3)]), color=c("blue","red")[mc$classification],
                                  main="Plot of the PC1, PC2 and PC3", pch=16)
     legend("right", legend=levels(as.factor(mc$classification)), col=c("blue", "red"), pch=16)
     screeplot(pca, npcs=20, type="lines", main="PCA screeplot")
     plot(with(pca,rotation)[,1:2],pch=19,cex=0.6)
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
  }
  pca_clustering_plot(pca,mc,paste0("~/Caprion/analysis/output/pca_clustering",suffix,".html"))
  list(pca=pca,km=km,mc=mc)
}

normalise_lr <- function(d,batches)
{
  caprion_lr <- d
  covars <- c("sexPulse","agePulse",paste0("PC",1:20))
  proteins <- setdiff(names(d),c("FID","IID","batch",pcs,covars))
  if (batches!=1) covars <- c(covars,pcs)
  mod <- model.matrix(as.formula(paste0("~",paste(covars,collapse="+"))), data=d)
  z <- sapply(proteins,
              function(col)
              {
                if ((batches==2 & gsub("^X([0-9])","\\1",col) %in% ZYQ.na) |
                    (batches==3 & gsub("^X([0-9])","\\1",col) %in% UDP.na)) {
                  x <- NA
                }
                else {
                  y <- invnormal(d[[col]])
                  l <- lm(y~mod[,-1])
                  r <- y-predict(l,na.action=na.pass)
                  x <- invnormal(r)
                }
                caprion_lr[col] <- x
                x
              })
  names(caprion_lr[proteins]) <- gsub("^X([0-9])","\\1",names(caprion_lr[proteins]))
  write.table(subset(caprion_lr[c("FID","IID")],!is.na(FID)),
              file=paste0("~/Caprion/analysis/output/caprion",suffix,"-",batches,".id"),
              quote=FALSE,col.names=FALSE,row.names=FALSE)
  write.table(subset(caprion_lr[c("FID","IID",proteins)],!is.na(FID)),
              file=paste0("~/Caprion/analysis/output/caprion",suffix,"-",batches,".mpheno"),
              quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
  write.table(subset(caprion_lr[c("FID","IID",proteins)],!is.na(FID)),
              file=paste0("~/Caprion/analysis/output/caprion",suffix,"-",batches,".pheno"),
              quote=FALSE,row.names=FALSE,sep="\t")
}

quantro_sparsenetgls <- function(edata,batch,mod)
{
  suppressMessages(library(quantro))
  suppressMessages(library(doParallel))
  registerDoParallel(cores=10)
  edata_quantro <- quantro(edata,batch,B=10000)
  print(edata_quantro)
  png(paste0("~/Caprion/analysis/output/matboxplot",suffix,".png"),width=12,height=10,unit="in",res=300)
  par(cex=0.6,cex.lab=0.6)
  matboxplot(edata,batch,pch=19)
  dev.off()
  png(paste0("~/Caprion/analysis/output/matdensity",suffix,".png"),width=12,height=10,unit="in",res=300)
  matdensity(edata, batch, pch=19, xlab = " ", ylab = "density", ylim=c(0,2),
             main = "Protein levels", brewer.n = 8, brewer.name = "Dark2")
  legend('top', c("ZWK", "ZYQ", "UDP"), col = c(1, 2, 3), lty = 1, lwd = 3)
  dev.off()
  if (FALSE)
  {
    suppressMessages(library(sparsenetgls))
    y <- t(edata)
    sngls <- sparsenetgls(y,mod)
    save(sngls,file=paste0("~/Caprion/analysis/output/sparsenetgls-",batch,suffix,".rda"))
    plotsngls(sngls,ith_lambda=5)
  }
}

options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(dplyr))
suppressMessages(library(gap))
load("~/Caprion/pilot/es.rda")

suffix <- "_dr"
ZYQ.na <- c("BROX","CT027","GHRL","PSB6")
UDP.na <- c("BROX","NCF2","SEM7A")
if (suffix=="_dr")
{
  protein_all <- protein_dr_all
  ZYQ.na <- c(ZYQ.na,"A1AG2","CAH13","COF1","CZIB","GBRL2","GLTD2","HPSE","ITA6","K2C5","MK14","NCF2",
                     "NDKA","NOMO1","NQO2","PTPRM","RGMA","SE6L2","SHLB1","SYYC","TMED8","VATG1","VAV")
  UDP.na <- c(UDP.na,"1433S", "AT1B1", "HPSE", "MYLK", "NENF", "NOMO1", "NQO2", "VATG1", "VAV")
}
prot <- t(exprs(protein_all))
colnames(prot) <- gsub("_HUMAN","",colnames(prot))

pnames <- colnames(prot)
xnames <- gsub("(^[0-9])","X\\1",pnames)
dnames <- colnames(prot)[grepl("^[0-9]",pnames)]
pcs <- paste0("ppc",1:3)
PCS <- paste0("PC",1:20)

ids <- c("identifier","Affymetrix_gwasQC_bl","caprion_id")
dates <- c("attendanceDate","monthPulse","yearPulse")
covars <- c("sexPulse","agePulse","ethnicPulse","processDate_bl","processTime_bl","classification")
data <- read.csv("~/Caprion/INTERVALdata_15SEP2021.csv")
eigenvec <- read.delim("~/Caprion/pilot/data/merged_imputation.eigenvec")
pilotsMap <- read.csv("~/Caprion/pilotsMap_15SEP2021.csv")
OmicsMap <- read.csv("~/Caprion/INTERVAL_OmicsMap_20210915.csv")
pca_km_mc <- pca_clustering(prot)
attach(pca_km_mc)
grouping <- data.frame(caprion_id=names(with(mc,classification)),classification=with(mc,classification)) %>%
            left_join(data.frame(with(pca,x)[,1:3],caprion_id=rownames(with(pca,x)))) %>%
            rename(ppc1=PC1,ppc2=PC2,ppc3=PC3)
detach(pca_km_mc)
id_data_covars <- merge(pilotsMap,OmicsMap,by="identifier",all=TRUE) %>%
                  merge(data,by="identifier",all=TRUE) %>%
                  left_join(grouping,by="caprion_id")
pheno <- select(OmicsMap,identifier,Affymetrix_gwasQC_bl,caprion_id) %>%
         merge(eigenvec,by.x="Affymetrix_gwasQC_bl",by.y="X.FID",all.x=TRUE) %>%
         left_join(grouping) %>%
         left_join(data) %>%
         right_join(data.frame(prot,caprion_id=rownames(prot))) %>%
         mutate(batch=match(substr(caprion_id,1,3),c("ZWK","ZYQ","UDP")))
write.table(select(pheno,Affymetrix_gwasQC_bl),file=paste0("~/Caprion/analysis/output/caprion",suffix,".id"),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
caprion_pheno <- mutate(pheno,FID=Affymetrix_gwasQC_bl,
                        IID=Affymetrix_gwasQC_bl)[c("FID","IID",ids,"batch",dates,covars,pcs,PCS,xnames)] %>%
                 arrange(batch,FID,IID)
names(caprion_pheno) <- gsub("^X([0-9])","\\1",names(caprion_pheno))
write.table(caprion_pheno,file=paste0("~/Caprion/analysis/output/caprion",suffix,".pheno"),
            quote=FALSE,row.names=FALSE,sep="\t")
allvars <- c("FID","IID","sexPulse","agePulse","batch",pcs,PCS,xnames)
dat <- mutate(pheno,FID=Affymetrix_gwasQC_bl,IID=Affymetrix_gwasQC_bl)[allvars]
mod <- model.matrix(as.formula(paste0(c("~sexPulse","agePulse",pcs,PCS),collapse="+")),
                    filter(pheno,!is.na(Affymetrix_gwasQC_bl)))
mcol <- apply(pheno[colnames(mod)[-1]],2,is.na)
many <- !apply(mcol,1,any)
for(batches in 1:3)
{
  d <- filter(dat[many,],batch==batches)
  z <- d
  names(z) <- gsub("^X([0-9])","\\1",names(z))
  write.table(z,file=paste0("~/Caprion/analysis/output/caprion",suffix,"-",batches,".tsv"),
              quote=FALSE,row.names=FALSE,sep="\t")
  normalise_lr(d,batches)
}

suppressMessages(library(sva))
edata <- protein_all
rownames(edata) <- sub("_HUMAN","",rownames(edata))
edata <- edata[!rownames(edata)%in%union(ZYQ.na,UDP.na),]
batch <- as.factor(match(substr(rownames(Biobase::pData(edata)),1,3),c("ZWK","ZYQ","UDP")))
table(batch)
# 1. parametric adjustment
combat_edata1 <- ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
# 2. non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
edata <- edata[,many]
batch <- pheno$batch[many]
quantro_sparsenetgls(edata,batch,mod)
# 3. reference-batch version, with covariates
combat_edata3 <- ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3, prior.plots=TRUE)
if (FALSE)
{
  pca_km_mc <- pca_clustering(combat_edata3)
  attach(pca_km_mc)
  if (interactive()) rgl::plot3d(with(pca,x[,c(2,1,3)]),col=mc$classification)
  pca_clustering_plot(pca,mc,paste0("~/Caprion/analysis/output/pca_clustering_combat",suffix,".html"))
  detach(pca_km_mc)
}
png(paste0("~/Caprion/analysis/output/combat-matboxplot",suffix,".png"),width=12,height=10,unit="in",res=300)
par(cex=0.6,cex.lab=0.6)
matboxplot(combat_edata3,groupFactor=batch, pch=19, ylab="Protein measurement")
dev.off()
png(paste0("~/Caprion/analysis/output/combat-matdensity",suffix,".png"),width=12,height=10,unit="in",res=300)
matdensity(combat_edata3,groupFactor=batch,xlab = " ", ylab = "density", pch=19, ylim=c(0,3),
           main = "Protein measurement", brewer.n = 8, brewer.name = "Dark2")
dev.off()
combat_quantro <- quantro(combat_edata3,batch,B=10000)
print(combat_quantro)
