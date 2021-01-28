#!/usr/bin/bash

export dir=~/rds/projects/olink_proteomics
R --no-save <<END
  library(openxlsx)
  dir <- Sys.getenv("dir")
# workbook
  wb <- file.path(dir,"ZYQ_EDR_28AUG2020.xlsx")
  Legend <- read.xlsx(wb, sheet = 1, startRow = 3)
  Samples <- read.xlsx(wb, sheet = 2, startRow = 6)
  names(Samples) <- c("LIMS.ID","Sample.ID","Comment")
  Mapping <- read.xlsx(wb, sheet = 3, startRow = 6)
  Annotations <- read.xlsx(wb, sheet = 4, startRow = 6)
# CSVs
  Comp_Neq1 <- read.csv(file.path(dir,"ZYQ_Comp_Neq1_Norm_Int_20200812.csv"))
  Normalized_All <- read.csv(file.path(dir,"ZYQ_Protein_Norm_All_20200813_v1.csv"))
  Protein_DR_Filt <- read.csv(file.path(dir,"ZYQ_Protein_Norm_DR_filt_20200813_v1.csv"))
  save(Legend,Samples,Mapping,Annotations,Comp_Neq1,Normalized_All,Protein_DR_Filt,file="2020.rda")
  extract2 <- function(prots=c("EPCR_HUMAN","PROC_HUMAN"))
  {
  # Piptides by Isotope.Group.ID
    mapping <- subset(Mapping,Protein%in%prots)
    d1 <- merge(mapping,Comp_Neq1,by="Isotope.Group.ID")
    piptides <- t(d1[grepl("Protein|ZYQ",names(d1))][,-1])
    colnames(piptides) <- paste0(gsub("HUMAN","",d1[,3]),d1[,1],"_",d1[,2])
  # Protein_All
    all <- Normalized_All[names(Normalized_All)%in%prots]
    colnames(all) <- gsub("HUMAN","All",colnames(all))
  # Protein_DR_Filt
    dr <- Protein_DR_Filt[names(Protein_DR_Filt)%in%prots]
    colnames(dr) <- gsub("HUMAN","DR",colnames(dr))
    cbind(piptides,all,dr)
  }
  load("2020.rda")
# EPCR-PROC
  epcr_proc <- extract2()
  library(corrplot)
  png(file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC-phase2-all.png"),width=12,height=10,units="in",pointsize=4,res=300)
  corrplot(cor(epcr_proc),method="square",type="lower",na.label="x")
  dev.off()
  library(pheatmap)
  short.mm <- c("PROC_442607670","PROC_442611348","PROC_442616032","PROC_442686905","PROC_442739582","PROC_442847259")
  r <- names(epcr_proc)
  r2 <- data.frame(long.name=r,short.name=substr(r,1,14))
  long.mm <- subset(r2,short.name%in%short.mm)[["long.name"]]
  png(file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC-phase2-corr.png"),width=12,height=10,units="in",pointsize=4,res=300)
  pheatmap(cor(epcr_proc[setdiff(names(epcr_proc),long.mm)]))
  dev.off()
#
  png(file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC-phase2-0.png"),width=12,height=10,units="in",pointsize=4,res=300)
  par(mfrow=c(3,1))
  with(epcr_proc,
  {
    plot(EPCR_All)
    hist(EPCR_All)
    boxplot(EPCR_All,horizontal = TRUE)
  })
  dev.off()
  png(file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC-phase2-1.png"),width=12,height=10,units="in",pointsize=4,res=300)
  par(mfrow=c(3,1))
  with(epcr_proc,
  {
    plot(PROC_All)
    hist(PROC_All)
    boxplot(PROC_All,horizontal = TRUE)
  })
  dev.off()
  epcr_proc_names <- colnames(epcr_proc)
  epcr_names <- epcr_proc_names[grepl("EPCR",epcr_proc_names)]
  EPCR <- data.frame(epcr_proc[,epcr_names])
  data.frame(names(EPCR))
  format(cor(EPCR),digits=3)
  head(EPCR)
  names(EPCR)[1:4] <- c("EPCR_442581804","EPCR_442582461","EPCR_442603139","EPCR_442605396")
  EPCR_lm <- lm(EPCR_All~EPCR_442581804+EPCR_442582461+EPCR_442603139+EPCR_442605396,data=EPCR)
  summary(EPCR_lm)

  require(ANN2)
  AE <- autoencoder(Normalized_All[,-1], hidden.layers=c(100,20,30), loss.type = 'pseudo-huber',
                    activ.functions = c('tanh','linear','tanh'),
                    batch.size = 8, optim.type = 'adam',
                    n.epochs = 1000, val.prop = 0)
# Plot loss during training
  plot(AE)
# Make reconstruction and compression plots
  reconstruction_plot(AE, Normalized_All[,-1])
  compression_plot(AE, Normalized_All[,-1])
# Reconstruct data and show states with highest anomaly scores
  recX <- reconstruct(AE, Normalized_All[,-1])
  sort(recX$anomaly_scores, decreasing = TRUE)
END
