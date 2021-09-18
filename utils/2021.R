# pilot, code ZWK

initialize <- function()
{
  rm(list=ls())
  library(Biobase)
  library(openxlsx)
  array_data <- function(data_frame,id,id_end_col)
  {
    rownames(data_frame) <- data_frame[[id]]
    arrayData <- as.matrix(data_frame[,-(1:id_end_col)])
  }
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
  peptides_ZWK <- ExpressionSet(peptideData,phenoData)
  save(protein_ZWK,dr_ZWK,peptides_ZWK,file="ZWK.rda")
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
  drData <- t(array_data(Protein_DR_Filt,"LIMS.ID",1))
  peptideData <- array_data(Comp_Neq1,"Isotope.Group.ID",5)
  all(rownames(Samples)==colnames(proteinData))
  protein_ZYQ <- ExpressionSet(proteinData,phenoData)
  dr_ZYQ <- ExpressionSet(drData,phenoData)
  peptides_ZYQ <- ExpressionSet(peptideData,phenoData)
  save(Legend,Samples,Mapping,Annotations,Comp_Neq1,Normalized_All,Protein_DR_Filt,protein_ZYQ,dr_ZYQ,peptides_ZYQ,file="ZYQ.rda")
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
}
