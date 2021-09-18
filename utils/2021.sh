#!/usr/bin/bash

export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot
if [ ! -d ${caprion}/data3 ]; then mkdir ${caprion}/data3; fi
if [ ! -d ${caprion}/bgen3 ]; then mkdir ${caprion}/bgen3; fi

R --no-save -q < ${caprion}/utils/2021.R

# UDP
R --no-save <<END
  caprion <- Sys.getenv("caprion")
  load("2021.rda")
  array_transpose <- function (x)
  {
    d <- x[,-1]
    rownames(d) <- x[,1]
    td <- t(d)
  }
  norm_all <- array_transpose(Protein_All_Peptides)
  dr_filt <- array_transpose(Protein_DR_Filt_Peptides)
  ppc <- prcomp(na.omit(norm_all), rank=10, scale=TRUE)
  pc1pc2 <- with(ppc,x)[,1:2]
  rownames(pc1pc2) <- rownames(norm_all)
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
  with(mc,
  {
       png(file.path(caprion,"data3","SERPING1.png"),res=300,width=12,height=10,units="in")
       scatterplot3d(with(ppc,x[,c(2,1,3)]), color=c("blue","red")[classification], main="Plot of the PC1, PC2 and PC3", pch=16)
       legend("right", legend=levels(as.factor(classification)), col=c("blue", "red"), pch=16)
       dev.off()
       plot3d(with(ppc,x[,c(2,1,3)]),col=classification)
  })
  pilotsMap <- read.csv("pilotsMap_15SEP2021.csv")
  OmicsMap <- read.csv("INTERVAL_OmicsMap_20210915.csv")
  data <- read.csv("INTERVALdata_15SEP2021.csv")
  id <- c("identifier","Affymetrix_gwasQC_bl","caprion_id")
  date <- c("attendanceDate","sexPulse","monthPulse","yearPulse","agePulse")
  covars <- c("ethnicPulse","ht_bl","wt_bl","CRP_bl","TRANSF_bl","processDate_bl","processTime_bl","classification")
  grouping <- data.frame(caprion_id=names(with(mc,classification)),classification=with(mc,classification))
  id_date_covars <- merge(merge(data,merge(pilotsMap,OmicsMap,by="identifier",all=TRUE),by="identifier",all=TRUE),grouping,by="caprion_id")
  samples <- merge(id_date_covars,Samples,by.x="caprion_id",by.y="LIMS.ID",all.y=TRUE) %>%
             left_join(pca,by=c("caprion_id"="id"))
  rownames(samples) <- samples[,1]
  library(Biobase)
  pap <- as.matrix(Protein_All_Peptides[,-1])
  rownames(pap) <- Protein_All_Peptides[,1]
  id_extra <- setdiff(rownames(samples),colnames(pap))
  annotations <- Annotations[,-1]
  rownames(annotations) <- Annotations[,1]
  experimentData <- new("MIAME", name="Contact", lab="Caprion", contact="contact@caprion",
         title="INTERVAL pilot", abstract="phase 2 ExpressionSet", url="email",
         other=list(notes="Created from csv files"))
  phenoData <- new("AnnotatedDataFrame", data=subset(samples,rownames(samples) %in% colnames(pap)))
  all(rownames(phenoData)==colnames(pap))
  library(pQTLtools)
# help("ExpressionSet-class")
  norm_es <- make_ExpressionSet(pap,phenoData,experimentData=experimentData)
  library(gap)
  r <- sapply(1:length(featureNames(norm_es)),function(r) {
                                               norm_es_r <- norm_es[r,]
                                               fn <- paste0("invnormal(",sub("(^[0-9])","X\\1",featureNames(norm_es_r)),")")
                                               f <- paste(fn,"~ agePulse + sexPulse + classification")
                                               z <- lm(as.formula(f),data=norm_es_r, na.action=na.exclude)
                                               resid(z)
                                               })
  colnames(r) <- featureNames(norm_es)
  d <- data.frame(r)
  d <- d %>%
       mutate(caprion_id=rownames(r)) %>%
       left_join(pData(norm_es)[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,setdiff(names(d),c("Affymetrix_gwasQC_bl","caprion_id"))) %>%
       filter(!caprion_id %in%c("UDP0138","UDP0481"))
  names(d) <- c("FID","IID",featureNames(norm_es))
  write.table(d,file=file.path(caprion,"data3","UDP.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
  checks <- function()
  {
    dim(pilotsMap)
    dim(OmicsMap)
    dim(data)
    head(pilotsMap)
    head(OmicsMap)
    head(data)
    head(exprs(norm_es))
    head(pData(phenoData))
    head(featureNames(norm_es))
    head(sampleNames(norm_es))
    experimentData(norm_es)
    intersect(OmicsMap$caprion_id,sampleNames(norm_es))
    intersect(OmicsMap$Affymetrix_gwasQC_bl,pData(norm_es)$Affymetrix_gwasQC_bl)
  nrows <- length(featureNames(norm_es))
  ncols <- length(intersect(OmicsMap$caprion_id,sampleNames(norm_es)))
  r <- matrix(NA,nrows,ncols)
  }
END
