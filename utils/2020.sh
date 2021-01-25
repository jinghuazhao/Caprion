#!/usr/bin/bash

export caprion=~/rds/projects/olink_proteomics
R --no-save <<END
  library(openxlsx)
  caprion <- Sys.getenv("caprion")
# workbook
  wb <- file.path(caprion,"ZYQ_EDR_28AUG2020.xlsx")
  Legend <- read.xlsx(wb, sheet = 1, startRow = 3)
  Samples <- read.xlsx(wb, sheet = 2, startRow = 6)
  names(Samples) <- c("LIMS.ID","Sample.ID","Comment")
  Mapping <- read.xlsx(wb, sheet = 3, startRow = 6)
  Annotations <- read.xlsx(wb, sheet = 4, startRow = 6)
# CSVs
  Normalized_All_Peptides <- read.csv(file.path(caprion,"ZYQ_Protein_Norm_All_20200813_v1.csv"))
  Protein_DR_Filt_Peptides <- read.csv(file.path(caprion,"ZYQ_Protein_Norm_DR_filt_20200813_v1.csv"))
  save(Legend,Samples,Annotations,Normalized_All_Peptides,Protein_DR_Filt_Peptides,file="2020.rda")
# EPCR
  EPCR <- Normalized_All_Peptides[c("LIMS.ID","EPCR_HUMAN")]
  EPCR_filt <- Protein_DR_Filt_Peptides[c("LIMS.ID","EPCR_HUMAN")]
# plot(cbind(EPCR,EPCR_filt)[,c(2,4)])
  par(mfrow=c(3,1))
  with(EPCR,
  {
    plot(EPCR_HUMAN)
    hist(EPCR_HUMAN)
    boxplot(EPCR_HUMAN,horizontal = TRUE)
  })
END
