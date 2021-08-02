#!/usr/bin/bash

function setup()
{
export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot
if [ ! -d ${caprion}/data3 ]; then mkdir ${caprion}/data3; fi
if [ ! -d ${caprion}/bgen3 ]; then mkdir ${caprion}/bgen3; fi

R --no-save <<END
  library(openxlsx)
  samples <- read.xlsx("UDP_EDR_20210423_samples.xlsx", sheet = 1, startRow = 5)
  wb <- "UDP_EDR_20210423.xlsx"
  Samples <- read.xlsx(wb,sheet="Samples",startRow=5)
  names(Samples) <- c("LIMS.ID","Sample.ID","Comment")
  Annotations <- read.xlsx(wb,sheet="Annotations",startRow=1)
  Mapping <- read.xlsx(wb,sheet="Mapping",startRow=6)
  Normalized_Peptides <- read.xlsx(wb,sheet="Normalized Peptides",startRow=1)
  Protein_All_Peptides <- read.xlsx(wb,sheet="Protein_All_Peptides",startRow=1)
  Protein_DR_Filt_Peptides <- read.xlsx(wb,sheet="Protein_DR_Filt_Peptides",startRow=1)
  save(Samples,Annotations,Mapping,Normalized_Peptides,Protein_All_Peptides,Protein_DR_Filt_Peptides,file="2021.rda")
END
}

function es()
{
R --no-save <<END
  caprion <- Sys.getenv("caprion")
  load("2021.rda")

  array_data <- function (x)
  {
    d <- x[,-1]
    rownames(d) <- x[,1]
    td <- t(d)
  }

  norm_all <- array_data(Protein_All_Peptides)
  dr_filt <- array_data(Protein_DR_Filt_Peptides)

  ppc <- prcomp(na.omit(norm_all), rank=50, scale=TRUE)
  pc1pc2 <- with(ppc,x)[,1:2]
  rownames(pc1pc2) <- rownames(norm_all)
  eigenvec <- with(ppc,rotation)[,1:2]
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
  pilotsMap <- read.csv("pilotsMap_17FEB2021.csv")
  OmicsMap <- read.csv("INTERVAL_OmicsMap_20210217.csv")
  data <- read.csv("INTERVALdata_17FEB2021.csv")
  head(pilotsMap)
  dim(pilotsMap)
  head(OmicsMap)
  dim(OmicsMap)
  head(data)
  dim(data)
  id <- c("identifier","Affymetrix_gwasQC_bl","caprion_id")
  date <- c("attendanceDate","sexPulse","monthPulse","yearPulse","agePulse")
  covars <- c("ethnicPulse","ht_bl","wt_bl","CRP_bl","TRANSF_bl","processDate_bl","processTime_bl","classification")
  grouping <- data.frame(caprion_id=names(with(mc,classification)),classification=with(mc,classification))
  id_date_covars <- merge(merge(data,merge(pilotsMap,OmicsMap,by="identifier",all=TRUE),by="identifier",all=TRUE),grouping,by="caprion_id")
  samples <- merge(id_date_covars,Samples,by.x="caprion_id",by.y="LIMS.ID",all.y=TRUE)
  rownames(samples) <- samples[,1]

# duplicates
# mapping <- Mapping[,-3]
# rownames(mapping) <- Mapping[,3]
# samples <- Samples[,-1]
# rownames(samples) <- Samples[,1]
  annotations <- Annotations[,-1]
  rownames(annotations) <- Annotations[,1]
  library(Biobase)
  experimentData <- new("MIAME", name="Contact", lab="Caprion", contact="contact@caprion",
         title="INTERVAL pilot", abstract="phase 2 ExpressionSet", url="email",
         other=list(notes="Created from csv files"))
  phenoData <- new("AnnotatedDataFrame", data=samples)
  library(pQTLtools)
  norm_es <- make_ExpressionSet(norm_all,phenoData,experimentData=experimentData)
  featureNames(norm_es)[1:10]
  sampleNames(norm_es)[1:10]
  experimentData(norm_es)
  lm(wt_bl~ht_bl,data=norm_es)
  lm(X1433B_HUMAN~wt_bl,data=norm_es)
END
}
