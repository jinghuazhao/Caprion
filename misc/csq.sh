#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr

function csq()
{
   Rscript -e '
     options(width=200)
     suppressMessages(require(dplyr))
     suppressMessages(require(pQTLtools))
     suppressMessages(require(stringr))
   # Caprion proteins
     caprion_cvt <- read.csv("~/Caprion/analysis/work/caprion.cis.vs.trans") %>%
                    dplyr::mutate(chr=SNPChrom,pos=SNPPos,MarkerName=SNP)
     caprion_dr_cvt <- read.csv("~/Caprion/analysis/work/caprion_dr.cis.vs.trans") %>%
                       dplyr::mutate(chr=SNPChrom,pos=SNPPos,MarkerName=SNP)
   # Caprion peptides
     peptide_cvt <- read.csv("~/Caprion/analysis/reports/peptide.cis.vs.trans") %>%
                    dplyr::mutate(chr=SNPChrom,pos=SNPPos,MarkerName=SNP)
     caprion <- caprion_cvt %>%
                select(chr,pos,MarkerName,prot,Type)
     caprion_dr <- caprion_dr_cvt %>%
                   select(chr,pos,MarkerName,prot,Type)
     peptide <- peptide_cvt %>%
                select(chr,pos,MarkerName,prot,Type)
   # VEP output
     vep <- "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/bgen/vep"
     pattern <- paste( protein_altering_variants, collapse = "|")
     suppressMessages(require(GenomicRanges))
     INF <- "/rds/project/rds-zuZwCZMsS0w/olink_proteomics/scallop/INF"
     plink <- "/rds/user/jhz22/hpc-work/bin/plink"
   # CSQ
     CSQ <- function(cvt)
     {
       b <- list()
       for (i in unique(dplyr::pull(cvt,chr))) {
         m <- dplyr::filter(cvt, chr %in% i) %>%
              dplyr::mutate(rsid = gsub("chr", "", MarkerName)) %>%
              dplyr::select(-MarkerName)
         u <- read.delim(file.path(vep, paste0("chr", i, ".tab.gz"))) %>%
              dplyr::select(Chrom, Pos, X.Uploaded_variation, Consequence) %>%
              setNames(c("chr", "pos", "rsid", "csq")) %>%
              dplyr::mutate(chr=as.character(chr),
                            chr=if_else(chr=="X","23",chr),chr=as.numeric(chr))
         bfile <- file.path(INF, "INTERVAL", "per_chr", paste0("snpid", i))
         b[[i]] <- csq(m, u, pattern, ldops = list(bfile = bfile, plink = plink))
       }
       dplyr::bind_rows(b) %>%
       dplyr::filter(r2 >= 0.8) %>%
       dplyr::rename(gene = prot) %>%
       dplyr::mutate(seqnames = as.integer(seqnames), pos = as.integer(pos)) %>%
       dplyr::arrange(seqnames, pos) %>%
       dplyr::select(-c(ref.seqnames, ref.start, ref.end, seqnames, pos, start, end)) %>%
       group_by(gene,rsid) %>%
       summarise(ref.rsid.all=paste(ref.rsid,collapse=";"),
                 ref.pos.all=paste(ref.pos,collapse=";"),
                 ref.csq.all=paste(ref.csq,collapse=";"),
                 r2.all=paste(r2,collapse=";"))
     }
     caprion_csq <- CSQ(caprion)
     save(caprion_cvt,caprion_csq,file="~/Caprion/analysis/reports/caprion_csq.rda",compress="xz",version=2)
     caprion_dr_csq <- CSQ(caprion_dr)
     save(caprion_dr_cvt,caprion_dr_csq,file="~/Caprion/analysis/reports/caprion_dr_csq.rda",compress="xz",version=2)
     peptide_csq <- CSQ(peptide)
     save(peptide_cvt,peptide_csq,file="~/Caprion/analysis/reports/peptide_csq.rda",compress="xz",version=2)
   # Examination
     library(psych)
     benchmarks <- c("A1BG","APOB","EPCR","ERAP2","PROC")
     load("~/Caprion/pilot/ZWK.rda")
     s <- data.frame()
     for (prot in benchmarks)
     {
       mps <- subset(mapping_ZWK,grepl(prot,Protein)) %>%
              mutate(mps=gsub("\\[.*?\\]", "-", Modified.Peptide.Sequence),width=nchar(mps))
       d <- describe(pull(mps,width)) %>%
            select(n,median,min,max,range,skew,kurtosis)
       rownames(d) <- prot
       s <- rbind(s,d)
     }
     knitr::kable(s,caption="Length of peptides",digits=2)
     load("~/Caprion/analysis/reports/caprion_csq.rda")
     dim(caprion)
     dim(caprion_csq)
     filter(caprion,prot%in%benchmarks)
     filter(caprion_csq,gene%in%benchmarks)
     caprion_csq_all <- left_join(caprion_csq,caprion,by=c('gene'='prot','rsid'='MarkerName'))
     table(caprion$Type)
     table(caprion_csq_all$Type)
     load("~/Caprion/analysis/reports/caprion_dr_csq.rda")
     dim(caprion_dr)
     dim(caprion_dr_csq)
     filter(caprion_dr,prot%in%benchmarks)
     filter(caprion_dr_csq,gene%in%benchmarks)
     caprion_dr_csq_all <- left_join(caprion_dr_csq,caprion_dr,by=c('gene'='prot','rsid'='MarkerName'))
     table(caprion_dr$Type)
     table(caprion_dr_csq_all$Type)
     load("~/Caprion/analysis/reports/peptide_csq.rda")
     dim(peptide)
     dim(peptide_csq)
     filter(peptide,prot%in%benchmarks)
     filter(peptide_csq,gene%in%benchmarks)
     peptide_csq_all <- left_join(peptide_csq,peptide,by=c('gene'='prot','rsid'='MarkerName'))
     table(peptide$Type)
     table(peptide_csq_all$Type)
  '
}
