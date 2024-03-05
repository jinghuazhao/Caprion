# ExpressionSet

rm(list=ls())
caprion <- Sys.getenv("caprion")
caprion <- ifelse(caprion=="",".",caprion)

suppressMessages(library(Biobase))
suppressMessages(library(openxlsx))
suppressMessages(library(pQTLtools))

array_data <- function(data_frame,id,id_end_col)
{
  arrayData <- data_frame[,-(1:id_end_col)]
  chars <- sapply(arrayData, is.character)
  arrayData[, chars] <- apply(arrayData[, chars], 2, as.numeric)
  rownames(arrayData) <- data_frame[[id]]
  as.matrix(arrayData)
}

zwk <- function()
# pilot, code ZWK
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
  mapping_ZWK <- rawIGs
  save(protein_ZWK,dr_ZWK,peptide_ZWK,mapping_ZWK,file="ZWK.rda")
}

zyq <- function()
# batch2, code ZYQ
{
  dir <- "/home/jhz22/Caprion/pre_qc_data/batch2/CAM1184-ZYQ"
  wb <- file.path(dir,"ZYQ_EDR_28AUG2020_Updated.xlsx")
  Legend <- read.xlsx(wb, sheet = 1, startRow = 3)
  Samples <- read.xlsx(wb, sheet = 2, startRow = 6)
  names(Samples) <- c("LIMS.ID","Sample.ID","Comment")
  Mapping <- read.xlsx(wb, sheet = 3, startRow = 6)
  Annotations <- read.xlsx(wb, sheet = 4, startRow = 6)
  Comp_Neq1 <- read.csv(file.path(dir,"ZYQ_Comp_Neq1_Norm_Int_20200812.csv"))
  Normalized_All <- read.csv(file.path(dir,"ZYQ_Protein_Norm_All_20200813_v1.csv"),check.names=FALSE)
  Protein_DR_Filt <- read.csv(file.path(dir,"ZYQ_Protein_Norm_DR_filt_20200813_v1.csv"),check.names=FALSE)
  rownames(Samples) <- Samples$"LIMS.ID"
  phenoData <- new("AnnotatedDataFrame", data=Samples)
  proteinData <- t(array_data(Normalized_All,"LIMS ID",1))
  drData <- t(array_data(Protein_DR_Filt,"LIMS ID",1))
  peptideData <- array_data(Comp_Neq1,"Isotope.Group.ID",5)
  all(rownames(Samples)==colnames(proteinData))
  protein_ZYQ <- ExpressionSet(proteinData,phenoData)
  dr_ZYQ <- ExpressionSet(drData,phenoData)
  peptide_ZYQ <- ExpressionSet(peptideData,phenoData)
  mapping_ZYQ <- Mapping
  save(protein_ZYQ,dr_ZYQ,peptide_ZYQ,mapping_ZYQ,file="ZYQ.rda")
}

udp <- function()
# batch3, code UDP
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
# help("ExpressionSet-class")
  library(pQTLtools)
  protein_UDP <- make_ExpressionSet(proteinData,phenoData,experimentData=experimentData)
  dr_UDP <- make_ExpressionSet(drData,phenoData,experimentData=experimentData)
  peptide_UDP <- make_ExpressionSet(peptideData,phenoData,experimentData=experimentData)
  mapping_UDP <- Mapping
  save(protein_UDP,dr_UDP,peptide_UDP,mapping_UDP,file="UDP.rda")
}

uhz <- function()
{
  dir <- "~/Caprion/pre_qc_data/batch4"
  wb <- file.path(dir,"UHZ_EDR_20240214_v1_Final.xlsx")
  Legend <- read.xlsx(wb, sheet=1, startRow=3)
  Samples <- read.xlsx(wb, sheet=2, startRow=5) %>%
             setNames(c("LIMS.ID","Sample.ID","Comment")) %>%
             mutate(LIMS.ID=gsub("L$","",LIMS.ID)) %>%
             `rownames<-`(.$LIMS.ID) %>%
             arrange(LIMS.ID)
  LIMS.ID <- pull(Samples,LIMS.ID)
  Annotations <- read.xlsx(wb, sheet=3, startRow=1)
  Mapping <- read.xlsx(wb, sheet=4, startRow=6)
  Normalized_Peptides <- read.xlsx(wb, sheet=5, startRow=1)
  Normalized_All <- read.xlsx(wb, sheet=6, startRow=1)
  Protein_DR_Filt <- read.xlsx(wb, sheet=7, startRow=1)
  Comp_Raw_Int <- read.csv(file.path(dir,"UHZ_Comp_Raw_Int_20240118_v1.csv"))
  phenoData <- AnnotatedDataFrame(Samples)
  proteinData <- array_data(Normalized_All,"Protein",1)
  drData <- array_data(Protein_DR_Filt,"Protein",1)
  peptideData <- array_data(Comp_Raw_Int[c("Isotope.Group.ID",LIMS.ID)],"Isotope.Group.ID",1)
  all(rownames(Samples)==colnames(proteinData))
  protein_UHZ <- ExpressionSet(proteinData,phenoData)
  dr_UHZ <- ExpressionSet(drData,phenoData)
  peptide_UHZ <- ExpressionSet(peptideData,phenoData)
  mapping_UHZ <- Mapping %>%
                 left_join(peptideData[1:5])
  save(protein_UHZ,dr_UHZ,peptide_UHZ,mapping_UHZ,file="UHZ.rda")
}

zwk();
zyq();
udp();
uhz();

contrast <- function()
{
  load("ZWK.rda")
  load("ZYQ.rda")
  load("UDP.rda")
  load("UHZ.rda")
  options(width=200)
  subset(mapping_ZWK[1:6],grepl("1433B",Protein)) %>% select(-Protein)
  subset(mapping_ZYQ,grepl("1433B",Protein)) %>% select(-Protein)
  subset(mapping_UDP,grepl("1433B",Protein)) %>% select(-Protein)
  subset(mapping_UHZ,grepl("1433B",Protein)) %>% select(-Protein)
}
