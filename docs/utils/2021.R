# pilot, code ZWK

rm(list=ls())
caprion <- Sys.getenv("caprion")

library(Biobase)
library(openxlsx)

array_data <- function(data_frame,id,id_end_col)
{
  rownames(data_frame) <- data_frame[[id]]
  arrayData <- as.matrix(data_frame[,-(1:id_end_col)])
}

zwk <- function()
{
  load("caprion.rda")
  rownames(Samples) <- Samples$"LIMS.ID:.Caprion.Sample.ID"
  phenoData <- new("AnnotatedDataFrame", data=Samples)
  proteinData <- array_data(Protein_All_Peptides,"Protein",1)
  drData <- array_data(Protein_DR_Filt_Peptides,"Protein",1)
  peptideData <- array_data(Normalized_Peptides,"Isotope.Group.ID",6)
  protein_ZWK <- ExpressionSet(proteinData,phenoData)
  dr_ZWK <- ExpressionSet(drData,phenoData)
  peptide_ZWK <- ExpressionSet(peptideData,phenoData)
  save(protein_ZWK,dr_ZWK,peptide_ZWK,file="ZWK.rda")
}

# batch2, code ZYQ

zyq <- function()
{
  dir <- "/home/jhz22/Caprion/pre_qc_data/batch2/CAM1184-ZYQ"
  wb <- file.path(dir,"ZYQ_EDR_28AUG2020_Updated.xlsx")
  Legend <- read.xlsx(wb, sheet = 1, startRow = 3)
  Samples <- read.xlsx(wb, sheet = 2, startRow = 6)
  names(Samples) <- c("LIMS.ID","Sample.ID","Comment")
  Mapping <- read.xlsx(wb, sheet = 3, startRow = 6)
  Annotations <- read.xlsx(wb, sheet = 4, startRow = 6)
  Comp_Neq1 <- read.csv(file.path(dir,"ZYQ_Comp_Neq1_Norm_Int_20200812.csv"))
  Normalized_All <- read.csv(file.path(dir,"ZYQ_Protein_Norm_All_20200813_v1.csv"))
  Protein_DR_Filt <- read.csv(file.path(dir,"ZYQ_Protein_Norm_DR_filt_20200813_v1.csv"))

  rownames(Samples) <- Samples$"LIMS.ID"
  phenoData <- new("AnnotatedDataFrame", data=Samples)
  proteinData <- t(array_data(Normalized_All,"LIMS.ID",1))
  rownames(proteinData) <- sub("^X","\\1",rownames(proteinData))
  drData <- t(array_data(Protein_DR_Filt,"LIMS.ID",1))
  rownames(drData) <- sub("^X","\\1",rownames(drData))
  peptideData <- array_data(Comp_Neq1,"Isotope.Group.ID",5)
  all(rownames(Samples)==colnames(proteinData))
  protein_ZYQ <- ExpressionSet(proteinData,phenoData)
  dr_ZYQ <- ExpressionSet(drData,phenoData)
  peptide_ZYQ <- ExpressionSet(peptideData,phenoData)
  save(Legend,Samples,Mapping,Annotations,Comp_Neq1,Normalized_All,Protein_DR_Filt,protein_ZYQ,dr_ZYQ,peptide_ZYQ,file="ZYQ.rda")
}

# batch3, code UDP
udp <- function()
{
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
# duplicates
# mapping <- Mapping[,-3]
# rownames(mapping) <- Mapping[,3]
# samples <- Samples[,-1]
# rownames(samples) <- Samples[,1]
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
  proteinData <- array_data(Protein_All_Peptides,"Protein",1)
  drData <- array_data(Protein_DR_Filt_Peptides,"Protein",1)
  peptideData <- array_data(Normalized_Peptides,"Isotope.Group.ID",6)
  id_extra <- setdiff(rownames(samples),colnames(proteinData))
  annotations <- Annotations[,-1]
  rownames(annotations) <- Annotations[,1]
  experimentData <- new("MIAME", name="Contact", lab="Caprion", contact="contact@caprion",
         title="INTERVAL pilot", abstract="batch3 ExpressionSet", url="email",
         other=list(notes="Created from csv files"))
  phenoData <- new("AnnotatedDataFrame", data=subset(samples,rownames(samples) %in% colnames(proteinData)))
  all(rownames(phenoData)==colnames(proteinData))
  library(pQTLtools)
# help("ExpressionSet-class")
  protein_UDP <- make_ExpressionSet(proteinData,phenoData,experimentData=experimentData)
  dr_UDP <- make_ExpressionSet(drData,phenoData,experimentData=experimentData)
  peptide_UDP <- make_ExpressionSet(peptideData,phenoData,experimentData=experimentData)
  save(protein_UDP,dr_UDP,peptide_UDP,file="UDP.rda")
  library(gap)
  r <- sapply(1:length(featureNames(protein_UDP)),function(r) {
                                               protein_UDP_r <- protein_UDP[r,]
                                               fn <- paste0("invnormal(",sub("(^[0-9])","X\\1",featureNames(protein_UDP_r)),")")
                                               f <- paste(fn,"~ agePulse + sexPulse + classification")
                                               z <- lm(as.formula(f),data=protein_UDP_r, na.action=na.exclude)
                                               resid(z)
                                               })
  colnames(r) <- featureNames(protein_UDP)
  d <- data.frame(r)
  d <- d %>%
       mutate(caprion_id=rownames(r)) %>%
       left_join(pData(protein_UDP)[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,setdiff(names(d),c("Affymetrix_gwasQC_bl","caprion_id"))) %>%
       filter(!caprion_id %in%c("UDP0138","UDP0481"))
  names(d) <- c("FID","IID",featureNames(protein_UDP))
  write.table(d,file=file.path(caprion,"data3","UDP.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
  save(protein_UDP,dr_UDP,peptide_UDP,file="UDP.rda")
  checks <- function()
  {
    dim(pilotsMap)
    dim(OmicsMap)
    dim(data)
    head(pilotsMap)
    head(OmicsMap)
    head(data)
    head(exprs(protein_UDP))
    head(pData(phenoData))
    head(featureNames(protein_UDP))
    head(sampleNames(protein_UDP))
    experimentData(protein_UDP)
    intersect(OmicsMap$caprion_id,sampleNames(protein_UDP))
    intersect(OmicsMap$Affymetrix_gwasQC_bl,pData(protein_UDP)$Affymetrix_gwasQC_bl)
  nrows <- length(featureNames(protein_UDP))
  ncols <- length(intersect(OmicsMap$caprion_id,sampleNames(protein_UDP)))
  r <- matrix(NA,nrows,ncols)
  }
}

overlap <- function(A,B) unlist(lapply(calculate.overlap(list(featureNames(A),featureNames(B))),length))

zwk();
zyq();
udp();
load("ZWK.rda")
load("ZYQ.rda")
load("UDP.rda")

library(VennDiagram)
overlap(protein_ZWK,protein_ZYQ)
overlap(protein_ZWK,protein_UDP)
overlap(protein_ZYQ,protein_UDP)
overlap(peptide_ZWK,peptide_ZYQ)
overlap(peptide_ZWK,peptide_UDP)
overlap(peptide_ZYQ,peptide_UDP)

protein <- list(pilot=featureNames(protein_ZWK),batch2=featureNames(protein_ZYQ),batch3=featureNames(protein_UDP))
dr <- list(pilot=featureNames(dr_ZWK),batch2=featureNames(dr_ZYQ),batch3=featureNames(dr_UDP))
peptide <- list(pilot=featureNames(peptide_ZWK),batch2=featureNames(peptide_ZYQ),batch3=featureNames(peptide_UDP))
unlist(lapply(calculate.overlap(protein),length))
unlist(lapply(calculate.overlap(dr),length))
unlist(lapply(calculate.overlap(peptide),length))
VennDiagram::venn.diagram(protein,"protein_ZWK-ZQY-UDP.png")
VennDiagram::venn.diagram(dr,"dr_ZWK-ZQY-UDP.png")
VennDiagram::venn.diagram(peptide,"peptide_ZWK-ZQY-UDP.png")
