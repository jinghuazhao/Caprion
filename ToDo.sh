#!/usr/bin/bash

pandoc ToDo.md --mathjax -s -o ToDo.html
st

export caprion=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics

export pre_qc_data=${caprion}/pre_qc_data
export peptide_ZWK=${pre_qc_data}/pilot
export peptide_ZYQ=${pre_qc_data}/batch2/CAM1184-ZYQ/
export peptide_UDP=${pre_qc_data}/batch3/CAM1184-UDP/

export ZYQ_R1=ZYQ_R1_Comp_Neq1_Raw_Int_Updated_20200505.txt
export ZYQ_R2=ZYQ_R2_Comp_Neq1_Raw_Int_Updated_20200606.txt
export ZYQ_R3=ZYQ_R3_Comp_Neq1_Raw_Int_Updated_20200608.txt
export ZYQ_R4=ZYQ_R4_Comp_Neq1_Raw_Int_Updated_20200615.txt

export UDP_R1=UDP_R1_Comp_Neq1_Raw_Int_Clean_20210412.txt
export UDP_R2=UDP_R2_Comp_Neq1_Raw_Int_Clean_20210412.txt
export UDP_R3=UDP_R3_Comp_Neq1_Raw_Int_Clean_20210412.txt

function header
{
  head -1 ${1} | \
  sed 's/Isotope Group ID/Isotope.Group.ID/;s/Intensity \[//g;s/\ (other)] (Pre-Combined Data)//g;s/\ (study)] (Pre-Combined Data)//g' > ${2}
}

header ${peptide_ZYQ}/${ZYQ_R1} ${caprion}/analysis/work/ZYQ_R1.header
header ${peptide_ZYQ}/${ZYQ_R2} ${caprion}/analysis/work/ZYQ_R2.header
header ${peptide_ZYQ}/${ZYQ_R3} ${caprion}/analysis/work/ZYQ_R3.header
header ${peptide_ZYQ}/${ZYQ_R4} ${caprion}/analysis/work/ZYQ_R4.header
header ${peptide_UDP}/${UDP_R1} ${caprion}/analysis/work/UDP_R1.header
header ${peptide_UDP}/${UDP_R2} ${caprion}/analysis/work/UDP_R2.header
header ${peptide_UDP}/${UDP_R3} ${caprion}/analysis/work/UDP_R3.header

Rscript -e '
  suppressMessages(library(dplyr))
  caprion <- Sys.getenv("caprion")
  peptide_ZYQ <- Sys.getenv("peptide_ZYQ")
  ZYQ_R1 <- read.delim(file.path(peptide_ZYQ,Sys.getenv("ZYQ_R1"))) %>%
            setNames(scan(file.path(caprion,"analysis","work","ZYQ_R1.header"),what="",sep="\t")) %>% mutate(R=1)
  ZYQ_R2 <- read.delim(file.path(peptide_ZYQ,Sys.getenv("ZYQ_R2"))) %>%
            setNames(scan(file.path(caprion,"analysis","work","ZYQ_R2.header"),what="",sep="\t")) %>% mutate(R=2)
  ZYQ_R3 <- read.delim(file.path(peptide_ZYQ,Sys.getenv("ZYQ_R3"))) %>%
            setNames(scan(file.path(caprion,"analysis","work","ZYQ_R3.header"),what="",sep="\t")) %>% mutate(R=3)
  ZYQ_R4 <- read.delim(file.path(peptide_ZYQ,Sys.getenv("ZYQ_R4"))) %>%
            setNames(scan(file.path(caprion,"analysis","work","ZYQ_R4.header"),what="",sep="\t")) %>% mutate(R=4)
  ZYQ_R <- bind_rows(ZYQ_R1,ZYQ_R2,ZYQ_R3,ZYQ_R4)
  save(ZYQ_R,file=file.path(caprion,"pilot","ZYQ_R.rda"))
  peptide_UDP <- Sys.getenv("peptide_UDP")
  UDP_R1 <- read.delim(file.path(peptide_UDP,Sys.getenv("UDP_R1"))) %>%
            setNames(scan(file.path(caprion,"analysis","work","UDP_R1.header"),what="",sep="\t")) %>% mutate(R=1)
  UDP_R2 <- read.delim(file.path(peptide_UDP,Sys.getenv("UDP_R2"))) %>%
            setNames(scan(file.path(caprion,"analysis","work","UDP_R2.header"),what="",sep="\t")) %>% mutate(R=2)
  UDP_R3 <- read.delim(file.path(peptide_UDP,Sys.getenv("UDP_R3"))) %>%
            setNames(scan(file.path(caprion,"analysis","work","UDP_R3.header"),what="",sep="\t")) %>% mutate(R=3)
  UDP_R <- bind_rows(UDP_R1,UDP_R2,UDP_R3)
  save(UDP_R,file=file.path(caprion,"pilot","UDP_R.rda"))
'

# /rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/pre_qc_data/spectra/
# /rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/pre_qc_data/spectral_library_ZWK/

Rscript -e '
options(width=200)
suppressMessages(library("dplyr"))
suppressMessages(library(Biobase))
setwd("../pilot")

mapping <- function(code,genes=c("PROC","EPCR","ERAP2"))
{
  if (code=="ZWK")
  {
     load("caprion.rda")
     dim(Protein_All_Peptides)
     dim(Protein_DR_Filt_Peptides)
     load("ZWK.rda")
     r <- data.frame(rawIGs,R=1)
     m <- Normalized_Peptides[c("Isotope.Group.ID","Modified.Peptide.Sequence","Protein")]
     prot <- protein_ZWK
     pept <- peptide_ZWK
  } else if (code=="ZYQ")
  {
     load("2020.rda")
     setdiff(colnames(Normalized_All),colnames(Protein_DR_Filt))
     dim(Normalized_All)
     dim(Protein_DR_Filt)
     load("ZYQ_R.rda")
     load("ZYQ.rda")
     r <- ZYQ_R
     m <- Mapping
     prot <- protein_ZYQ
     pept <- peptide_ZYQ
  } else
  {
     load("2021.rda")
     dim(Protein_All_Peptides)
     dim(Protein_DR_Filt_Peptides)
     setdiff(Protein_All_Peptides$Protein,Protein_DR_Filt_Peptides$Protein)
     dim(Normalized_Peptides)
     load("UDP_R.rda")
     load("UDP.rda")
     r <- UDP_R
     m <- Mapping
     prot <- protein_UDP
     pept <- peptide_UDP
  }
  z <- subset(m,grepl("PROC|EPCR|ERAP2",Protein)) %>%
       arrange(desc(Protein),Modified.Peptide.Sequence)
  print(table(z$Protein))
  igi_PROC <- subset(m,grepl("PROC",Protein)) %>%
              pull(Isotope.Group.ID)
  igi_EPCR <- subset(m,grepl("EPCR",Protein)) %>%
              pull(Isotope.Group.ID)
  igi_ERAP2 <- subset(m,grepl("ERAP2",Protein)) %>%
               pull(Isotope.Group.ID)
  igi <- c(igi_PROC,igi_EPCR,igi_ERAP2)
  r_n <- names(r)
  r_ZWK <- r_n[grepl("ZWK",r_n)]
  print(r_ZWK)
  x <- filter(r,Isotope.Group.ID %in% igi) %>%
       select(Isotope.Group.ID,r_n[grepl(code,r_n)],R)
  rd <- aggregate(select(x,-c(Isotope.Group.ID,R)),by=list(x$Isotope.Group.ID),max,na.rm=TRUE) %>%
         rename(Isotope.Group.ID=Group.1)
  rd[rd==-Inf] <- NA
  write.csv(rd,file=paste0("~/",code,".csv"),quote=FALSE,row.names=FALSE)
  diff_PROC <- setdiff(igi_PROC,intersect(rd$Isotope.Group.ID,igi_PROC))
  diff_EPCR <- setdiff(igi_EPCR,intersect(rd$Isotope.Group.ID,igi_EPCR))
  diff_ERAP2 <- setdiff(igi_ERAP2,intersect(rd$Isotope.Group.ID,igi_ERAP2))
  diff_all <- setdiff(igi,intersect(rd$Isotope.Group.ID,igi))
  cat(code,"diff_PROC:",diff_PROC,"\n")
  cat(code,"diff_EPCR:",diff_EPCR,"\n")
  cat(code,"diff_ERAP2:",diff_ERAP2,"\n")
  cat(code,"diff_All:",diff_all,"\n")
  dr <- apply(rd>50000,1,sum)/ncol(rd)*100
  d <- data.frame(Isotope.Group.ID=as.numeric(rd$Isotope.Group.ID),dr=dr) %>%
       left_join(m) %>%
       arrange(Protein)
  igi_pept <- intersect(paste(igi),featureNames(pept))
  prot_pept <- combine(prot[c("PROC_HUMAN","EPCR_HUMAN","ERAP2_HUMAN")],pept[igi_pept])
  r_a <- cor(t(rd[-1]),use="everything")
  colnames(r_a) <- rownames(r_a) <- d$Isotope.Group.ID
  rd_b <- rd
  rd_b[rd_b==0] <- NA
  r_b <- cor(t(rd_b[-1]),use="everything")
  colnames(r_b) <- rownames(r_b) <- d$Isotope.Group.ID
  rd_c <- rd
  rd_c[!is.na(rd_c)&rd_c<50000] <- NA
  r_c <- cor(t(rd_c[-1]),use="everything")
  colnames(r_c) <- rownames(r_c) <- d$Isotope.Group.ID
  rd_d <- t(exprs(pept[igi_pept]))
  r_d <- cor(rd_d,use="everything")
  colnames(r_d) <- rownames(r_d) <- d$Isotope.Group.ID
  rd_e <- t(exprs(pept[igi_pept]))
  rd_e[rd_e==0] <- NA
  r_e <- cor(rd_e,use="everything")
  rd_f <- t(exprs(pept[igi_pept]))
  rd_f[!is.na(t(rd_c[-1]))&t(rd_c[-1])<50000] <- NA
  r_f <- cor(rd_f,use="everything")
  r <- cor(t(exprs(prot_pept)),use="everything")
  colnames(r) <- rownames(r) <- c("PROC_HUMAN","EPCR_HUMAN","ERAP2_HUMAN",d$Isotope.Group.ID)
  invisible(list(z=z,d=rd,diff_PROC=diff_PROC,diff_EPCR=diff_EPCR,diff_ERAP2=diff_ERAP2,diff_all=diff_all,dr=d,
                 r_a=r_a,r_b=r_b,r_c=r_c,r_d=r_d,r_e=r_e,r_f=r_f,r=r))
}

zwk <- mapping("ZWK")
zyq <- mapping("ZYQ")
udp <- mapping("UDP")

dr2 <- function(a,b,c)
# difference in r^2
{
  all3 <- intersect(intersect(colnames(a),colnames(b)),colnames(c))
  ab <- sum(abs(a[all3,all3]^2-b[all3,all3]^2),na.rm=TRUE)
  ac <- sum(abs(a[all3,all3]^2-c[all3,all3]^2),na.rm=TRUE)
  bc <- sum(abs(b[all3,all3]^2-c[all3,all3]^2),na.rm=TRUE)
  c(ab,ac,bc)
}

dr2(zwk$r_a,zyq$r_a,udp$r_a)
dr2(zwk$r_b,zyq$r_b,udp$r_b)
dr2(zwk$r_c,zyq$r_c,udp$r_c)
dr2(zwk$r_d,zyq$r_d,udp$r_d)
dr2(zwk$r_e,zyq$r_e,udp$r_e)
dr2(zwk$r_f,zyq$r_f,udp$r_f)
dr2(zwk$r,zyq$r,udp$r)

wr2 <- function(a,b,c)
# weight in r^2
{
  all3 <- intersect(intersect(colnames(a),colnames(b)),colnames(c))
  N <- length(all3)
  wa <- sum(a[all3,all3]^2,na.rm=TRUE) - N
  wb <- sum(b[all3,all3]^2,na.rm=TRUE) - N
  wc <- sum(c[all3,all3]^2,na.rm=TRUE) - N
  c(N,wa,wb,wc)
}

wr2(zwk$r_a,zyq$r_a,udp$r_a)
wr2(zwk$r_b,zyq$r_b,udp$r_b)
wr2(zwk$r_c,zyq$r_c,udp$r_c)
wr2(zwk$r_d,zyq$r_d,udp$r_d)
wr2(zwk$r_e,zyq$r_e,udp$r_e)
wr2(zwk$r_f,zyq$r_f,udp$r_f)
wr2(zwk$r,zyq$r,udp$r)
'
