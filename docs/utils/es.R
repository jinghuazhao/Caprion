load("2020.rda")

array_data <- function (x)
{
  d <- x[,-1]
  rownames(d) <- x[,1]
  td <- t(d)
  rownames(td) <- gsub("^X","",rownames(td))
  print(dim(td))
  td
}

norm_all <- array_data(Normalized_All)
dr_filt <- array_data(Protein_DR_Filt)

ppc <- with(Normalized_All, prcomp(na.omit(Normalized_All[,-1]), rank=50, scale=TRUE))
pc1pc2 <- with(ppc,x)[,1:2]
rownames(pc1pc2) <- Normalized_All[["LIMS.ID"]]
eigenvec <- with(ppc,rotation)[,1:2]
library(mclust)
mc <- Mclust(pc1pc2,G=2)
summary(mc)
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
