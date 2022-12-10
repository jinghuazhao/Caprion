#!/usr/bin/bash

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

function peptides()
{
  cd ~/Caprion/pilot/bgen2/EPCR-PROC
  cat <(echo Protein Group Seq OBS_CT BETA SE T_STAT P | tr ' ' '\t') \
      <(zgrep -w rs867186 PROC* EPCR* | \
        sed 's/_invn-plink2.gz//;s/:/\t/' | \
        sed 's/_/\t/g;s/All/All\t/;s/DR/DR\t/' | \
        sort -k1,1 -k3,3 | \
        cut -f4-10 --complement) > ll
}

# /rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/pre_qc_data/spectra/
# /rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/pre_qc_data/spectral_library_ZWK/

Rscript -e '
options(width=200)
suppressMessages(library("dplyr"))
suppressMessages(library(Biobase))
suppressMessages(library(pheatmap))
caprion <- Sys.getenv("caprion")
setwd("../pilot")

mapping <- function(code,genes=c("PROC","EPCR","ERAP2"))
{
  if (code=="ZWK")
  {
     load("caprion.rda")
     dim(Protein_All_Peptides)
     dim(Protein_DR_Filt_Peptides)
     load("ZWK.rda")
     r <- data.frame(rawIGs,R=1) %>%
          select(-Protein)
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
  r_other <- setdiff(r_n[!grepl(paste("ZWK",code,sep="|"),r_n)][-(1:5)],"R")
  print(r_other)
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
       arrange(Protein) %>%
       mutate(Protein=gsub("_HUMAN","",Protein))
  igi_pept <- intersect(paste(igi),featureNames(pept))
  igi_other <- filter(r,Isotope.Group.ID%in%igi_pept)[1:5]
  names(igi_other) <- paste0(code,"_",names(igi_other))
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
  rd_f[!is.na(t(rd[-1]))&t(rd[-1])<50000] <- NA
  r_f <- cor(rd_f,use="everything")
  r_pp <- cor(t(exprs(prot_pept)),use="everything")
  colnames(r_pp) <- rownames(r_pp) <- c("PROC_HUMAN","EPCR_HUMAN","ERAP2_HUMAN",d$Isotope.Group.ID)
  m <- subset(z,Isotope.Group.ID%in%d$Isotope.Group.ID) %>%
      mutate(Protein=gsub("_HUMAN","",Protein),code=code)
  d <- mutate(d,code=code)
  invisible(list(z=z,r=r,d=rd,m=m,pp=prot_pept,igi_other=unique(igi_other),
                 diff_PROC=diff_PROC,diff_EPCR=diff_EPCR,diff_ERAP2=diff_ERAP2,diff_all=diff_all,dr=d,
                 r_a=r_a,r_b=r_b,r_c=r_c,r_d=r_d,r_e=r_e,r_f=r_f,r_pp=r_pp))
}

zwk <- mapping("ZWK")
zyq <- mapping("ZYQ")
udp <- mapping("UDP")

q1 <- bind_rows(zwk$m,zyq$m,udp$m) %>%
      group_by(Protein,Modified.Peptide.Sequence,code) %>%
      summarise(groups=paste(Isotope.Group.ID,collapse=";")) %>%
      mutate(protein_peptide_isotope=paste(Protein,Modified.Peptide.Sequence,groups,sep="_"))
t1  <- with(q1,table(protein_peptide_isotope,code))
knitr::kable(t1,caption="Protein/Peptide/Isotope_Group by batch")
write.csv(t1,file=file.path("~/Q1.csv"),quote=FALSE)

q2 <- bind_rows(zwk$dr,zyq$dr,udp$dr) %>%
      filter(dr>=10) %>%
      group_by(Protein,Modified.Peptide.Sequence,code) %>%
      summarise(groups=paste(Isotope.Group.ID,collapse=";")) %>%
      mutate(protein_peptide_isotope=paste(Protein,Modified.Peptide.Sequence,groups,sep="_"))
t2 <- with(q2,table(protein_peptide_isotope,code))
knitr::kable(t2,caption="Protein/Peptide/Isotope_Group by batch")
write.csv(t2,file=file.path("~/Q2.csv"),quote=FALSE)

q3 <- function(dlist,src="peptide",rt="EPCR-PROC-ERAP2-corr")
{
  m <- with(dlist,m)
  d <- with(dlist,d)
  rownames(d) <- d[[1]]
  cpmi <- select(m,code,Protein,Modified.Peptide.Sequence,Isotope.Group.ID) %>%
          mutate(name=paste(code,Protein,Modified.Peptide.Sequence,Isotope.Group.ID,sep="_"))
  if (src=="peptide")
  {
    fn <- setdiff(featureNames(pp),c("EPCR_HUMAN","ERAP2_HUMAN","PROC_HUMAN"))
    ann <- left_join(data.frame(Isotope.Group.ID=as.numeric(fn)),cpmi) %>%
           arrange(Protein,Modified.Peptide.Sequence) %>%
           select(Isotope.Group.ID,name)
    dat <- pp[paste(ann$Isotope.Group.ID)]
    featureNames(dat) <- ann$name
    corr <- cor(t(exprs(dat)),use="everything")
  } else if (src=="intensity") {
    ann <- left_join(d[1],cpmi) %>%
           arrange(Protein,Modified.Peptide.Sequence) %>%
           select(Isotope.Group.ID,name)
    dat <- t(d[-1][paste(ann$Isotope.Group.ID),])
    colnames(dat) <- ann$name
    corr <- cor(dat,use="everything")
  }
  f_png <- file.path(caprion,"analysis","work",paste0(rt,".png"))
  png(f_png,width=12,height=10,units="in",pointsize=4,res=300)
  pheatmap(corr,cluster_rows=FALSE, cluster_cols=FALSE)
  dev.off()
  corr[!lower.tri(cor(corr))] <- NA
  f_csv <- file.path(caprion,"analysis","work",paste0(rt,".csv"))
  write.csv(format(corr,digits=2),file=f_csv,quote=FALSE)
}

q3(zwk,src="intensity",rt="intensity_ZWK-corr")
q3(zyq,src="intensity",rt="intensity_ZYQ-corr")
q3(udp,src="intensity",rt="intensity_UDP-corr")

q3(zwk,src="peptide",rt="PP_ZWK-corr")
q3(zyq,src="peptide",rt="PP_ZYQ-corr")
q3(udp,src="peptide",rt="PP_UDP-corr")

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
dr2(zwk$r_pp,zyq$r_pp,udp$r_pp)

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
wr2(zwk$r_pp,zyq$r_pp,udp$r_pp)

other_check <- left_join(zwk$igi_other,zyq$igi_other,by=c('ZWK_Isotope.Group.ID'='ZYQ_Isotope.Group.ID')) %>%
               left_join(udp$igi_other,by=c('ZWK_Isotope.Group.ID'='UDP_Isotope.Group.ID')) %>%
               rename(Isotope.Group.ID=ZWK_Isotope.Group.ID)
other_check_n <- names(other_check)
n1 <- other_check_n[grepl("Isotope.Group.ID|Modified.Peptide.Sequence",other_check_n)]
n2 <- other_check_n[grepl("Monoisotopic",other_check_n)]
n3 <- other_check_n[grepl("Time",other_check_n)]
n4 <- other_check_n[grepl("Charge",other_check_n)]
knitr::kable(other_check[c(n1,n3,n4,n2)])
)

knitr::kable(other_check[c("Isotope.Group.ID",other_check_n[grepl("Modified.Peptide.Sequence",other_check_n)])])
knitr::kable(other_check[c("Isotope.Group.ID",other_check_n[grepl("Monoisotopic",other_check_n)])])
knitr::kable(other_check[c("Isotope.Group.ID",other_check_n[grepl("Time",other_check_n)])])
knitr::kable(other_check[c("Isotope.Group.ID",other_check_n[grepl("Charge",other_check_n)])])
'

#!/usr/bin/bash

#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --job-name=miamiplot
#SBATCH --array=317,751
#SBATCH --qos=cardio
#SBATCH --mem=50000
#SBATCH --time=12:00:00
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/work/_miamiplot_%A_%a.o
#SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/work/_miamiplot_%A_%a.e

export TMPDIR=${HPC_WORK}/work
export caprion=~/Caprion/pilot
export analysis=~/Caprion/analysis

function miamiplot2()
{
  export TMPDIR=/rds/user/jhz22/hpc-work/work
  export caprion=~/rds/projects/Caprion_proteomics/pilot
  export analysis=~/Caprion/analysis
  export prot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${caprion}/2020.id | sed -r 's/^X([0-9])/\1/')
  export uniprot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' ${caprion}/2020.id)
  module load gcc/6
  if [ ! -f ${analysis}/work/${uniprot}-${prot}-all.dat.gz ]; then
     (
       echo snpid chr pos rsid p z
       gunzip -c ${caprion}/bgen2/EPCR-PROC/${prot}_All_invn-plink2.gz | \
       awk 'NR>1{if($4<$5) {a1=$4;a2=$5} else {a1=$5;a2=$4}; snpid=$1":"$3"_"a1"_"a2; print snpid,$1,$2,$3,$12,$11}' | \
       sort -k2,2n -k3,3n
     ) | \
     gzip -f > ${analysis}/work/${uniprot}-${prot}-all.dat.gz
  fi
  if [ ! -f ${analysis}/work/${uniprot}-${prot}-dr.dat.gz ]; then
     (
       echo snpid chr pos rsid p z
       gunzip -c ${caprion}/bgen2/EPCR-PROC/${prot}_DR_invn-plink2.gz  | \
       awk 'NR>1{if($4<$5) {a1=$4;a2=$5} else {a1=$5;a2=$4}; snpid=$1":"$3"_"a1"_"a2; print snpid,$1,$2,$3,$12,$11}' | \
       sort -k2,2n -k3,3n
     ) | \
     gzip -f > ${analysis}/work/${uniprot}-${prot}-dr.dat.gz
  fi
  Rscript -e '
  caprion <- Sys.getenv("caprion");analysis <- Sys.getenv("analysis"); prot <- Sys.getenv("prot"); uniprot <- Sys.getenv("uniprot")
  protein <- prot
  suppressMessages(require(dplyr))
  suppressMessages(require(gap))
  gwas1 <- read.table(file.path(analysis,"work",paste(uniprot,prot,"all.dat.gz",sep="-")),as.is=TRUE,header=TRUE) %>% filter(!is.na(p))
  gwas2 <- read.table(file.path(analysis,"work",paste(uniprot,prot,"dr.dat.gz",sep="-")),as.is=TRUE,header=TRUE) %>% filter(!is.na(p))
  png(file.path(analysis,"work",paste(uniprot,prot,"all-dr.png",sep="-")),res=300,width=12,height=10,units="in")
  chrmaxpos <- miamiplot2(gwas1,gwas2,name1="All",name2="DR",z1="z",z2="z") %>%
               filter(chr!=-Inf & maxpos!=-Inf & genomestartpos!=-Inf & labpos!=-Inf)
  labelManhattan(chr=20,pos=33764554,name=rs867186,gwas1,gwasZLab="z",chrmaxpos=chrmaxpos)
  dev.off()
  '
# rm ${analysis}/work/${uniprot}-${prot}-all.dat.gz
# rm ${analysis}/work/${uniprot}-${prot}-dr.dat.gz
}

miamiplot2

# ---

for batch in ZWK ZYQ UDP
do
  for src in intensity PP
  do
    convert -resize 50% work/${src}_${batch}-corr.png ${src}_${batch}-corr.png
  done
done
convert -resize 50% work/Q9UNN8-EPCR-all-dr.png EPCR-All-DR.png
convert -resize 50% work/P04070-PROC-all-dr.png PROC-All-DR.png
pandoc ToDo.md --mathjax -s -o ToDo.html
st
