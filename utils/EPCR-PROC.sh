#!/usr/bin/bash

export caprion=~/rds/projects/Caprion_proteomics/pilot/
export protein=EPCR-PROC
export snp_pos=~/rds/projects/olink_proteomics/scallop/INF/work/snp_pos

(
  echo MarkerName Chrom Start End P prot rsid
  echo ---------- ----- ----- --- - ---- ----
  grep Q9UNN8 ${caprion}/1e-5/caprion-invn.sentinels | sort -k6,6 | join -16 -22 - ${snp_pos} | \
  awk '{gsub("chr","",$2)};1' | \
  sort -k2,2n -k3,3n
) | cut -d ' ' -f2-4 --complement | sed 's/ /|/g' > ${caprion}/EPCR-PROC/Q9UNN8_invn.sentinels
(
  echo MarkerName Chrom Start End P prot rsid
  echo ---------- ----- ----- --- - ---- ----
  grep P04070 ${caprion}/1e-5/caprion-invn.sentinels | sort -k6,6 | join -16 -22 - ${snp_pos} | \
  awk '{gsub("chr","",$2)};1' | \
  sort -k2,2n -k3,3n
) | cut -d' ' -f2-4 --complement | sed 's/ /|/g' > ${caprion}/EPCR-PROC/P04070_invn.sentinels

qpdf ${caprion}/scatter-histogram-boxwhisker.pdf --pages . 319,754 -- ${caprion}/${protein}/${protein}-desc.pdf
convert ${caprion}/${protein}/${protein}-desc.pdf ${caprion}/${protein}/${protein}-desc.png

# net use W: \\ME-FILER1\GROUPS$\MGEU
# W:\Factors\INTERVAL\Caprion_proteomics

R --no-save <<END
options(width=500)
library(pQTLtools)

annotation <- function()
## Annotations
{
  library(pQTLtools)
  pag <- subset(caprion,Accession %in% c("Q9UNN8","P04070"))
#    Protein Accession  Gene            Protein.Description
# EPCR_HUMAN    Q9UNN8 PROCR Endothelial protein C receptor
# PROC_HUMAN    P04070  PROC Vitamin K-dependent protein C

  subset(hg19Tables,hgncSym %in% c("PROCR","PROC"))
# chrom chromStart chromEnd strand    acc uniprotName geneName geneSynonyms hgncSym         ensGene
# chr20   33759957 33764613      + Q9UNN8  EPCR_HUMAN    PROCR         EPCR   PROCR ENSG00000101000
#  chr2 128177518 128186519      + P04070  PROC_HUMAN    PROC                  PROC ENSG00000115718

  qqman <- pag[c("Accession","Gene")]
  write.table(qqman,file="qqman.list",quote=FALSE,row.names=FALSE,col.names=FALSE)
}

load("caprion.rda")

extract <- function(prots=c("EPCR_HUMAN","PROC_HUMAN"))
{
# Isotope.Group.ID, Modified.Peptide.Sequence
  d1 <- subset(Normalized_Peptides,Protein %in% prots)
  piptides <- t(d1[,-(1:6)])
  colnames(piptides) <- paste0(gsub("HUMAN","",d1[,2]),d1[,1],"_",d1[,3])
# Protein
  d2 <- subset(Protein_All_Peptides,Protein %in% prots)
  all <- t(d2[,-1])
  colnames(all) <- gsub("HUMAN","All",d2[,1])
# Protein_DR_Filt
  d3 <- subset(Protein_All_Peptides,Protein %in% prots)
  dr <- t(d3[,-1])
  colnames(dr) <- gsub("HUMAN","DR",d3[,1])
  cbind(piptides,all,dr)
}

epcr_proc <- extract()
dir <- Sys.getenv("caprion")
write.table(epcr_proc,file=file.path(dir,"EPCR-PROC","EPCR-PROC.tsv"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
EPCR_PROC_corr <- cor(epcr_proc)
library(pheatmap)
png(file.path(dir,"EPCR-PROC","EPCR-PROC-corr.png"),width=12,height=10,units="in",pointsize=4,res=300)
pheatmap(cor(epcr_proc))
dev.off()
EPCR_PROC_corr[!lower.tri(cor(epcr_proc))] <- NA
write.table(format(EPCR_PROC_corr,digits=2),file=file.path(dir,"EPCR-PROC","EPCR-PROC-corr.tsv"),quote=FALSE,sep="\t")

epcr_proc_names <- colnames(epcr_proc)
epcr_names <- epcr_proc_names[grepl("EPCR",epcr_proc_names)]
EPCR <- data.frame(epcr_proc[,epcr_names])
data.frame(names(EPCR))
format(cor(EPCR),digits=3)
names(EPCR)[1:4] <- c("EPCR_442603139","EPCR_442605396","EPCR_442582461","EPCR_442581804")
EPCR_lm <- lm(EPCR_All~EPCR_442603139+EPCR_442605396+EPCR_442582461+EPCR_442581804,data=EPCR)
summary(EPCR_lm)

END

cut -d' ' -f2,3 --complement --output-delimiter='|' \
${caprion}/bgen2/EPCR-PROC/5e-8/caprion-invn.sentinels | \
sed 's/Chrom/chr/;s/_invn//g;s/chr[0-9]*://' > ${caprion}/EPCR-PROC/EPCR-PROC.sentinels

R --no-save -q <<END
  caprion <- Sys.getenv("caprion")
  sentinels <- read.table(file.path(caprion,"EPCR-PROC","EPCR-PROC.sentinels"),header=TRUE,sep="|")
  ldmatrix <- TwoSampleMR::ld_matrix(with(sentinels,SNP),with_alleles=FALSE)^2
  ldmatrix[upper.tri(ldmatrix)] <- NA
  round(ldmatrix,digits=2)
END
