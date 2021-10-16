# Phenotypic information on UDP

caprion <- Sys.getenv("caprion")
caprion <- ifelse(caprion=="",".",caprion)
load("UDP.rda")

library(Biobase)
library(dplyr)
library(gap)

checks <- function()
{
  dim(pilotsMap)
  dim(OmicsMap)
  dim(data)
  head(pilotsMap)
  head(OmicsMap)
  head(data)
  head(exprs(protein_UDP))
  head(featureNames(protein_UDP))
  head(sampleNames(protein_UDP))
  experimentData(protein_UDP)
  intersect(OmicsMap$caprion_id,sampleNames(protein_UDP))
  intersect(OmicsMap$Affymetrix_gwasQC_bl,pData(protein_UDP)$Affymetrix_gwasQC_bl)
  nrows <- length(featureNames(protein_UDP))
  ncols <- length(intersect(OmicsMap$caprion_id,sampleNames(protein_UDP)))
  r <- matrix(NA,nrows,ncols)
}

cluster <- function(data, interactive=FALSE)
{
  ppc <- prcomp(na.omit(data), rank=10, scale=TRUE)
  pc1pc2 <- with(ppc,x)[,1:2]
  rownames(pc1pc2) <- rownames(data)
  eigenvec <- with(ppc,rotation)[,1:2]
  library(dplyr)
  pca <- with(ppc,x) %>%
         data.frame()
  pca <- pca %>%
         mutate(id=rownames(pca))
  library(mclust)
  mc <- Mclust(pc1pc2,G=2)
  summary(mc)
  library(scatterplot3d)
  library(rgl)
  if (interactive) with(mc,
  {
       png(file.path(caprion,"data3","3d.png"), res=300, width=12, height=10, units="in")
       scatterplot3d(with(ppc,x[,c(2,1,3)]), color=c("blue","red")[classification], main="Plot of the PC1, PC2 and PC3", pch=16)
       legend("right", legend=levels(as.factor(classification)), col=c("blue", "red"), pch=16)
       dev.off()
       plot3d(with(ppc,x[,c(2,1,3)]),col=classification)
  })
  mc
}

UDP <- function()
{
  pilotsMap <- read.csv("pilotsMap_15SEP2021.csv")
  OmicsMap <- read.csv("INTERVAL_OmicsMap_20210915.csv")
  data <- read.csv("INTERVALdata_15SEP2021.csv")
  norm_all <- t(exprs(protein_UDP))
  id <- c("identifier","Affymetrix_gwasQC_bl","caprion_id")
  date <- c("attendanceDate","sexPulse","agePulse")
  covars <- c("ethnicPulse","ht_bl","wt_bl","classification")
  mc <- cluster(norm_all)
  grouping <- data.frame(caprion_id=names(with(mc,classification)),classification=with(mc,classification))
  id_date_covars <- merge(merge(data,merge(pilotsMap,OmicsMap,by="identifier",all=TRUE),by="identifier",all=TRUE),grouping,by="caprion_id",all=TRUE)
  ev20 <- read.delim(file.path(caprion,"data","merged_imputation.eigenvec"))
  names(ev20)[1] <- "FID"
  samples <- merge(id_date_covars,data.frame(caprion_id=sampleNames(protein_UDP)),by="caprion_id",all.y=TRUE) %>%
             left_join(ev20,by=c("Affymetrix_gwasQC_bl"="FID")) %>%
             select(caprion_id,sexPulse,agePulse,ethnicPulse,ht_bl,wt_bl,Affymetrix_gwasQC_bl,classification,paste0("PC",1:20)) %>%
             left_join(data.frame(pData(protein_UDP)),by=c("caprion_id"="LIMS.ID..Caprion.Sample.ID"))
  rownames(samples) <- samples[["caprion_id"]]
  pheno <- data.frame(caprion_id=sampleNames(protein_UDP),t(exprs(protein_UDP))) %>%
           right_join(samples,by="caprion_id")
  rhs <- paste("agePulse","sexPulse","classification",paste0("PC",1:20,collapse="+"),sep="+")
  r <- sapply(featureNames(protein_UDP),function(r) {
              x <- sub("(^[0-9])","X\\1",r)
              y <- paste0("invnormal(",x,")")
              f <- paste(y,"~",rhs)
              z <- lm(formula(f),data=pheno[c(x,"agePulse","sexPulse","classification",paste0("PC",1:20))], na.action=na.exclude)
              resid(z)
            })
  rownames(r) <- sampleNames(protein_UDP)
  d <- data.frame(r) %>%
       mutate(caprion_id=rownames(r)) %>%
       left_join(pheno[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,names(data.frame(r))) %>%
       filter(!caprion_id %in%c("UDP0138","UDP0481")) %>%
       mutate(caprion_id=Affymetrix_gwasQC_bl)
  names(d) <- c("FID","IID",gsub("HUMAN","invn",featureNames(protein_UDP)))
  write.table(d,file=file.path(caprion,"data3","UDP.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
}

UDP()

# --- legacy code ---

array_transpose <- function (x)
{
  d <- x[,-1]
  rownames(d) <- x[,1]
  td <- t(d)
}
