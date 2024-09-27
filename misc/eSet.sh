#!/usr/bin/bash

#SBATCH --job-name=_utils
#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --error=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/work/eSet.e
#SBATCH --output=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/work/eSet.o

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
export PERL5LIB=
module load ceuadmin/R
module load samtools/1.13/gcc/zwxn7ug3
module load perl/5.26.3_system/gcc-8.4.1-4cl2czq
module load libiconv/1.16/intel/64iicvbf
module load ceuadmin/ensembl-vep/111-icelake

export TMPDIR=${HPC_WORK}/work
export analysis=~/Caprion/analysis
export pre_qc_data=/rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics
export suffix=_dr

function eSet()
{
  Rscript -e '
     options(width = 200)
     pre_qc_data <- Sys.getenv("pre_qc_data")
     analysis <- Sys.getenv("analysis")
     pilot <- gsub("analysis", "pilot", analysis)

     raw_process <- function(f)
     {
       n <- read.delim(f,nrows=1) %>% length
       h <- scan(f,sep="\t",what="",n=n)
       h <- gsub(".*\\[(UDP\\d+|XCC\\w+|ZWK\\w+|ZYQ\\w+).*\\].*", "\\1",h)
       read.table(f,col.names=h,sep="\t",skip=1)
     }
     require(openxlsx)
     suppressMessages(library(dplyr))
     ZWK_path <- file.path(pre_qc_data,"pilot","ZWK_EDR_20191002.xlsx")
     raw_ZWK <- openxlsx::read.xlsx(ZWK_path, sheet="Raw IGs",colNames=TRUE, skipEmptyRows=TRUE, startRow=5)
     neq_ZWK <- openxlsx::read.xlsx(ZWK_path, sheet="Normalized Peptides",colNames=TRUE, skipEmptyRows=TRUE, startRow=5)

     ZYQ_path <- file.path(pre_qc_data,"batch2","CAM1184-ZYQ")
     f1 <- file.path(ZYQ_path,"ZYQ_R1_Comp_Neq1_Raw_Int_Updated_20200505.txt")
     f2 <- file.path(ZYQ_path,"ZYQ_R3_Comp_Neq1_Raw_Int_Updated_20200608.txt")
     f3 <- file.path(ZYQ_path,"ZYQ_R2_Comp_Neq1_Raw_Int_Updated_20200606.txt")
     f4 <- file.path(ZYQ_path,"ZYQ_R4_Comp_Neq1_Raw_Int_Updated_20200615.txt")
     raw_ZYQ1 <- raw_process(f1)
     raw_ZYQ2 <- raw_process(f2)
     raw_ZYQ3 <- raw_process(f3)
     raw_ZYQ4 <- raw_process(f4)
     raw_ZYQ <- raw_ZYQ1[names(raw_ZYQ1)[grep("UDP|Isotope.Group.ID", names(raw_ZYQ1))]] %>%
                left_join(raw_ZYQ2[names(raw_ZYQ2)[grep("ZYQ|Isotope.Group.ID", names(raw_ZYQ2))]]) %>%
                left_join(raw_ZYQ3[names(raw_ZYQ3)[grep("ZYQ|Isotope.Group.ID", names(raw_ZYQ3))]]) %>%
                left_join(raw_ZYQ4[names(raw_ZYQ4)[grep("ZYQ|Isotope.Group.ID", names(raw_ZYQ4))]])
     neq_ZYQ <- read.csv(file.path(ZYQ_path,"ZYQ_Comp_Neq1_Norm_Int_20200812.csv"))

     UDP_path <- file.path(pre_qc_data,"batch3","CAM1184-UDP")
     f1 <- file.path(UDP_path,"UDP_R1_Comp_Neq1_Raw_Int_Clean_20210412.txt")
     f2 <- file.path(UDP_path,"UDP_R2_Comp_Neq1_Raw_Int_Clean_20210412.txt")
     f3 <- file.path(UDP_path,"UDP_R3_Comp_Neq1_Raw_Int_Clean_20210412.txt")
     raw_UDP1 <- raw_process(f1)
     raw_UDP2 <- raw_process(f2)
     raw_UDP3 <- raw_process(f3)
     raw_UDP <- raw_UDP1[names(raw_UDP1)[grep("UDP|Isotope.Group.ID", names(raw_UDP1))]] %>%
                left_join(raw_UDP2[names(raw_UDP2)[grep("UDP|Isotope.Group.ID", names(raw_UDP2))]]) %>%
                left_join(raw_UDP3[names(raw_UDP3)[grep("UDP|Isotope.Group.ID", names(raw_UDP3))]])
     neq_UDP <- openxlsx::read.xlsx(file.path(UDP_path,"UDP_EDR_20210423.xlsx"),
                                    sheet="Normalized Peptides",colNames=TRUE, skipEmptyRows=TRUE, startRow=1)
     save(list=ls(),file=file.path(analysis,"work","eSet.rda"))
  '
}

eSet
