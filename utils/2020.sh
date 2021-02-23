#!/usr/bin/bash

export dir=~/rds/projects/olink_proteomics
export caprion=~/rds/projects/olink_proteomics/scallop/Caprion
if [ ! -d ${dir}/scallop/Caprion/data2 ]; then mkdir ${dir}/scallop/Caprion/data2; fi
R --no-save <<END
  dir <- Sys.getenv("dir")
# workbook
  library(openxlsx)
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
    peptides_all_dr <- cbind(piptides,all,dr)
    write.table(data.frame(caprion_id=rownames(peptides_all_dr),peptides_all_dr),file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC_All.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    peptides_all_dr
  }
  load("2020.rda")
# EPCR-PROC
  epcr_proc <- extract2()
  library(corrplot)
  png(file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC-phase2-all.png"),width=12,height=10,units="in",pointsize=4,res=300)
  EPCR_PROC_corr <- cor(epcr_proc)
  corrplot(EPCR_PROC_corr,method="square",type="lower",na.label="x")
  dev.off()
  library(pheatmap)
  short.mm <- c("PROC_442607670","PROC_442611348","PROC_442616032","PROC_442686905","PROC_442739582","PROC_442847259")
  r <- names(epcr_proc)
  r2 <- data.frame(long.name=r,short.name=substr(r,1,14))
  long.mm <- subset(r2,short.name%in%short.mm)[["long.name"]]
  EPCR_PROC_corr <- cor(epcr_proc[setdiff(names(epcr_proc),long.mm)])
  png(file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC-phase2-corr.png"),width=12,height=10,units="in",pointsize=4,res=300)
  pheatmap(EPCR_PROC_corr)
  dev.off()
  EPCR_PROC_corr[!lower.tri(EPCR_PROC_corr)] <- NA
  write.table(format(EPCR_PROC_corr,digits=2),file=file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC-phase2-corr.tsv"),quote=FALSE,sep="\t")
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
# Make phenotype file
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
  covars <- c("ethnicPulse","ht_bl","wt_bl","CRP_bl","TRANSF_bl","processDate_bl","processTime_bl")
  id_date_covars <- merge(data,merge(pilotsMap,OmicsMap,by="identifier",all=TRUE),by="identifier",all=TRUE)
  dim(id_date_covars)
  head(id_date_covars[c(id,date,covars)])
  peptides_all_dr <- read.delim(file.path(dir,"scallop","Caprion","EPCR-PROC","EPCR-PROC_All.tsv"),as.is=TRUE)
  pheno2 <- merge(id_date_covars[c(id,date,covars)],peptides_all_dr,by="caprion_id")
  write.table(subset(pheno2,!is.na(Affymetrix_gwasQC_bl),select=Affymetrix_gwasQC_bl),
              file=file.path(dir,"scallop","Caprion","data2","affymetrix.id"),quote=FALSE,row.names=FALSE,col.names=FALSE)
  library(gap)
  C <- "agePulse"
  D <- "sexPulse"
  cols <- grep("EPCR|PROC",names(pheno2))
  P <- names(pheno2)[cols]
  P_inv <- sapply(cols,function(x) {invnormal(pheno2[,x])})
  colnames(P_inv) <- paste0(P,"_invn")
  pheno2 <- within(data.frame(pheno2,P_inv),{ID_1 <- Affymetrix_gwasQC_bl; ID_2 <- Affymetrix_gwasQC_bl; missing <- 0})
  snptest_sample(subset(pheno2,!is.na(Affymetrix_gwasQC_bl)),sample_file=file.path(dir,"scallop","Caprion","data2","caprion.sample"),
                 ID_1 = "ID_1",ID_2 = "ID_2", missing = "missing", C = C, D = D, P = paste0(P,"_invn"))
END
(
  cut -d' ' -f3-5 --complement ${caprion}/data2/caprion.sample | awk '{if(NR==1) {$1="FID";$2="IID"}};1' | sed '2d' > ${caprion}/data2/caprion.pheno
  cut -d' ' -f1,2,4,5 ${caprion}/data2/caprion.sample | awk '{if(NR==1) {$1="FID";$2="IID"}};1' | sed '1,2d' > ${caprion}/data2/caprion.covar
)
