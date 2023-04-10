#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export batch=1
export SLURM_ARRAY_TASK_ID=14
export protein=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${pilot}/work/caprion.varlist)
export interval=${HPC_WORK}/data/interval

Rscript -e '
 library(Biobase)
 library(dplyr)
 batches <- c("ZWK","ZYQ","UDP")
 batch <- Sys.getenv("batch") %>% as.numeric
 code <- batches[batch]
 protein <- Sys.getenv("protein")
 load(paste0("~/Caprion/pilot/",code,".rda"))
 peptide_ids <- subset(get(paste("mapping",code,sep="_")),Protein==paste0(protein,"_HUMAN"))[["Isotope.Group.ID"]]
 peptide_exprs <- exprs(get(paste("peptide",code,sep="_")))
 peptide_exprs <- subset(peptide_exprs, rownames(peptide_exprs) %in% peptide_ids)
 peptide_dat <- data.frame(id=colnames(peptide_exprs),t(peptide_exprs))
 write.table(peptide_dat,file=paste0("~/",code,".pheno"),col.names=gsub("^[X]","",names(peptide_dat)),row.names=FALSE,quote=FALSE)
'
