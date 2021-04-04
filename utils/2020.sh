#!/usr/bin/bash

export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot
if [ ! -d ${caprion}/data2 ]; then mkdir ${caprion}/data2; fi
if [ ! -d ${caprion}/bgen2 ]; then mkdir ${caprion}/bgen2; fi

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
  load("2020.rda")
  Annotations[,1] <- gsub("1433","X1433",Annotations[,1])
  Annotations[,1] <- gsub("4F","X4F",Annotations[,1])
  Annotations[,1] <- gsub("6P","X6P",Annotations[,1])
  prot_uniprot <- data.frame(prot=gsub("_HUMAN","",Annotations[,1]),Accession=Annotations[,2])
  write.table(prot_uniprot,file="2020.id",col.names=FALSE,quote=FALSE,row.names=FALSE)
# AutoEncoder
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
# PCA
  ppc <- with(Normalized_All, prcomp(na.omit(Normalized_All[,-1]), rank=50, scale=TRUE))
  pc1pc2 <- with(ppc,x)[,1:2]
  rownames(pc1pc2) <- Normalized_All[["LIMS.ID"]]
  eigenvec <- with(ppc,rotation)[,1:2]
  cor(eigenvec)
  cor(pc1pc2)
  pdf("pca-2020.pdf")
  screeplot(ppc, npcs=20, type="lines", main="PCA screeplot")
  plot(eigenvec,pch=19,cex=0.6)
  title("Eigenvectors")
  plot(pc1pc2,pch=19,cex=0.6)
  title("Principal components")
  biplot(ppc,cex=0.1)
  title("biplot")
  dev.off()
  pdf("clustering-2020.pdf")
# K-means clustering
  km <- kmeans(pc1pc2,2)
  table(with(km,cluster))
  plot(pc1pc2, col = with(km,cluster), pch=19, cex=0.8)
  points(with(km,centers), col = 1:2, pch = 8, cex = 2)
  title("K-means clustering")
# Model-based clustering
  library(mclust)
  mc <- Mclust(pc1pc2,G=2)
  summary(mc)
  table(with(mc,classification))
  plot(mc, what=c("classification"))
  title("Model-based clustering")
  dev.off()
  caprion_mc <- read.csv("ZYQ_PC1_groups_20200703.csv")
  mc_caprion_mc <- cbind(caprion_mc,classification=with(mc,classification))
  with(mc_caprion_mc,table(pc1_group,classification))
# Phenotype files
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
  dim(id_date_covars)
  head(id_date_covars[c(id,date,covars)])
# All data at phase 2
  extract <- function()
  {
  # Piptides by Isotope.Group.ID
    d1 <- merge(Mapping,Comp_Neq1,by="Isotope.Group.ID")
    piptides <- t(d1[grepl("Protein|ZYQ",names(d1))][,-1])
    colnames(piptides) <- paste0(gsub("HUMAN","",d1[,3]),d1[,1],"_",d1[,2])
    piptides <- data.frame(caprion_id=rownames(piptides),piptides)
  # Protein_All
    all <- Normalized_All[names(Normalized_All)]
    colnames(all) <- gsub("HUMAN","All",colnames(all))
  # Protein_DR_Filt
    dr <- Protein_DR_Filt[names(Protein_DR_Filt)]
    colnames(dr) <- gsub("HUMAN","DR",colnames(dr))
    All <- merge(merge(piptides,all,by.x="caprion_id",by.y="LIMS.ID"),dr,by.x="caprion_id",by.y="LIMS.ID")
    write.table(All,file.path(dir,"pilot","data2","phase2_All.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
    All
  }
# SNPTEST phenotype file
  missing <- read.table("data/merged_imputation.missing",col.names=c("affymetrix_gwasqc_bl","missing"))
  eigenvec <- read.delim("data/merged_imputation.eigenvec")
  missing_eigenvec <- merge(missing,eigenvec,by.x="affymetrix_gwasqc_bl",by.y="X.FID")
  library(gap)
  peptides_all_dr <- extract()
  id_date_covars_missing_eigenvec <- merge(within(id_date_covars[c(id,date,covars)],{bmi=wt_bl/ht_bl/ht_bl;ethnicity=1}),
                                           missing_eigenvec,
                                           by.x="Affymetrix_gwasQC_bl",by.y="affymetrix_gwasqc_bl",all.x=TRUE)
  nonwhite <- with(id_date_covars_missing_eigenvec,ethnicPulse%in%c("Not Disclosed","Unknown"))
  id_date_covars_missing_eigenvec[nonwhite,"ethnicity"] <- 2
  id_date_covars_missing_eigenvec_peptides_all_dr <- merge(id_date_covars_missing_eigenvec,peptides_all_dr,by="caprion_id")
  ord <- with(id_date_covars_missing_eigenvec_peptides_all_dr,order(Affymetrix_gwasQC_bl))
  pheno2 <- id_date_covars_missing_eigenvec_peptides_all_dr[ord,]
  id1_id2_missing <- with(pheno2, data.frame(ID_1=Affymetrix_gwasQC_bl, ID_2=Affymetrix_gwasQC_bl, missing=missing))
  C <- c("agePulse","bmi",paste0("PC",1:20))
  D <- c("sexPulse","ethnicity","classification")
  CD <- pheno2[c(C,D)]
  cols <- 41:ncol(pheno2)
  P <- names(pheno2)[cols]
  P_invn <- sapply(cols,function(x) {invnormal(pheno2[,x])})
  colnames(P_invn) <- paste0(P,"_invn")
  write.table(subset(pheno2,!is.na(Affymetrix_gwasQC_bl),select=Affymetrix_gwasQC_bl),
              file=file.path(dir,"pilot","data2","affymetrix.id"),quote=FALSE,row.names=FALSE,col.names=FALSE)
  snptest_sample(subset(data.frame(id1_id2_missing,CD,P_invn),!is.na(ID_1)),
                 sample_file=file.path(dir,"pilot","data2","phase2.sample"),
                 ID_1 = "ID_1",ID_2 = "ID_2", missing = "missing", C = C, D = D, P = paste0(P,"_invn"))
# EPCR-PROC
  subset(Annotations,grepl("EPCR",Protein)|grepl("PROC",Protein))
  vars <- names(peptides_all_dr)
  prots <- c("EPCR","PROC")
  epcr_proc <- peptides_all_dr[grepl("caprion_id",vars)|grepl(prots[1],vars)|grepl(prots[2],vars)]
  write.table(epcr_proc,file.path(dir,"pilot","EPCR-PROC","EPCR-PROC_All.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
  epcr_proc_names <- names(epcr_proc)[-1]
  epcr_names <- epcr_proc_names[grepl("EPCR",epcr_proc_names)]
  EPCR <- data.frame(epcr_proc[,epcr_names])
  data.frame(names(EPCR))
  format(cor(EPCR),digits=3)
  head(EPCR)
  names(EPCR) <- substr(epcr_names,1,14)
  options(width=100)
  round(cor(EPCR[c("EPCR_442581804","EPCR_442582461","EPCR_442603139","EPCR_442605396","EPCR_All","EPCR_DR")]),digits=3)
  EPCR_lm <- lm(EPCR_All~EPCR_442581804+EPCR_442582461+EPCR_442603139+EPCR_442605396,data=EPCR)
  proc_names <- epcr_proc_names[grepl("PROC",epcr_proc_names)]
  PROC <- data.frame(epcr_proc[,proc_names])
  data.frame(names(PROC))
  format(cor(PROC),digits=3)
  head(PROC)
  names(PROC) <- substr(proc_names,1,14)
  library(corrplot)
  png(file.path(dir,"pilot","EPCR-PROC","EPCR-PROC-phase2-all.png"),width=12,height=10,units="in",pointsize=4,res=300)
  EPCR_PROC_corr <- cor(epcr_proc[,-1])
  corrplot(EPCR_PROC_corr,method="square",type="lower",na.label="x")
  dev.off()
  library(pheatmap)
  png(file.path(dir,"pilot","EPCR-PROC","EPCR-PROC-phase2-corr.png"),width=12,height=10,units="in",pointsize=4,res=300)
  col_row <- rownames(subset(as.data.frame(EPCR_PROC_corr),!is.na(EPCR_All)))
  pheatmap(EPCR_PROC_corr[col_row,col_row])
  dev.off()
  write.table(format(EPCR_PROC_corr[col_row,col_row],digits=2),file=file.path(dir,"pilot","EPCR-PROC","EPCR-PROC-phase2-corr.tsv"),quote=FALSE,sep="\t")
  png(file.path(dir,"pilot","EPCR-PROC","EPCR-PROC-phase2-0.png"),width=12,height=10,units="in",pointsize=4,res=300)
  par(mfrow=c(3,1))
  with(epcr_proc,
  {
    plot(EPCR_All)
    hist(EPCR_All)
    boxplot(EPCR_All,horizontal = TRUE)
  })
  dev.off()
  png(file.path(dir,"pilot","EPCR-PROC","EPCR-PROC-phase2-1.png"),width=12,height=10,units="in",pointsize=4,res=300)
  par(mfrow=c(3,1))
  with(epcr_proc,
  {
    plot(PROC_All)
    hist(PROC_All)
    boxplot(PROC_All,horizontal = TRUE)
  })
  dev.off()
# SNPTEST phenotype file
  epcr_proc_P_invn <- P_invn[,colnames(P_invn)%in%paste0(epcr_proc_names,"_invn")]
  snptest_sample(subset(data.frame(id1_id2_missing,CD,epcr_proc_P_invn),!is.na(ID_1)),
                 sample_file=file.path(dir,"pilot","data2","epcr-proc.sample"),
                 ID_1 = "ID_1",ID_2 = "ID_2", missing = "missing", C = C, D = D, P = paste0(epcr_proc_names,"_invn"))
END

# All for phase 2
(
  cut -d' ' -f3-28 --complement ${caprion}/data2/phase2.sample | awk '
     {if(NR==1) {$1="FID";$2="IID"} else gsub(/NA/,"-999",$0)};1' | \
     sed '2d' > ${caprion}/data2/phase2.pheno
  cut -d' ' -f1,2,4-28 ${caprion}/data2/phase2.sample | awk '{if(NR==1) {$1="FID";$2="IID"}};1' | sed '1,2d' > ${caprion}/data2/phase2.covar
  sed '1,2d' ${caprion}/data2/phase2.sample | awk '$28==1 {print $1,$2}' > ${caprion}/data2/phase2.group1
  sed '1,2d' ${caprion}/data2/phase2.sample | awk '$28==2 {print $1,$2}' > ${caprion}/data2/phase2.group2
)
# EPCR-PROC
(
  cut -d' ' -f3-28 --complement ${caprion}/data2/epcr-proc.sample | awk '
      {if(NR==1) {$1="FID";$2="IID"} else gsub(/NA/,"-999",$0)};1' | \
      sed '2d' > ${caprion}/data2/epcr-proc.pheno
  cut -d' ' -f1,2,4-28 ${caprion}/data2/epcr-proc.sample | awk '{if(NR==1) {$1="FID";$2="IID"}};1' | sed '1,2d' > ${caprion}/data2/epcr-proc.covar
  sed '1,2d' ${caprion}/data2/epcr-proc.sample | awk '$28==1 {print $1,$2}' > ${caprion}/data2/epcr-proc.group1
  sed '1,2d' ${caprion}/data2/epcr-proc.sample | awk '$28==2 {print $1,$2}' > ${caprion}/data2/epcr-proc.group2
)

# Comp_Neq1
# Normalized log-intensities for every study sample (identified by their Caprion's LIMS ID), and identified IG. IGs mapping to several proteins were not normalized and are therefore not reported.
# An IG is a peak (from Rosetta Elucidator), identified by its Isotope Group ID, matched or not to a Modified Peptide Sequence, along with the corresponding Monoisotopic m/z, retention time (Max Isotope Time Centroid) and Charge

# Normalised log-intensities mapped to single proteins for All samples
export All=$(head ${caprion}/data2/phase2_All.tsv | sed 's/\t/\n/g' | grep "_All$")

# Normalised log-intensities mapped to single proteins for samples with Detection Rate > 10%
# DR is calculated as the proportion of samples with a raw intensity (i.e. pre-normalization, original - not logged - scale) above 50,000.
export DR=$(head ${caprion}/data2/phase2_All.tsv | sed 's/\t/\n/g' | grep "_DR$")

# a list of unplotted Miami plot
ls miamiplot/|sed 's/-phase1-phase2.png//;s/-/\t/' | cut -f2 | grep -f - -v 2020.id > 2020.left

export threshold=5e-8
export a=${caprion}/bgen/${threshold}/caprion-invn.sentinels
export b=${caprion}/bgen2/${threshold}/caprion-invn.sentinels
bedtools intersect -a <(awk '{if(NR>1) $1="chr"$1;print}' $a | tr ' ' '\t') \
                   -b <(awk '{if(NR>1) $1="chr"$1;print}' $b | tr ' ' '\t') \
                   -wa -wb -loj | \
awk '$8!="."'
