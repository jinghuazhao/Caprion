#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis

for SLURM_ARRAY_TASK_ID in 14
do
export SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}
export protein=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${pilot}/work/caprion.varlist)
if [ ! -d ${analysis}/peptide/${protein} ]; then mkdir ${analysis}/peptide/${protein}; fi
export sbatch=${analysis}/peptide/${protein}/${protein}.sb

Rscript -e '
  suppressMessages(library(Biobase))
  suppressMessages(library(dplyr))
  pilot <- Sys.getenv("pilot")
  analysis <- Sys.getenv("analysis")
  protein <- Sys.getenv("protein")
  pcs <- paste0("ppc",1:3)
  ccovars <- c("sexPulse","agePulse",pcs,paste0("PC",1:20))
  dat <- read.delim(file.path(pilot,"work","caprion.pheno"))[c("FID","IID","caprion_id",ccovars)] %>%
         rename(id=caprion_id)
  normalise_peptide <- function(batch,batches=c("ZWK","ZYQ","UDP"),verbose=FALSE)
  {
    if (batch==1) covars <- setdiff(ccovars,pcs) else covars <- ccovars
    code <- batches[batch]
    load(paste0("~/Caprion/pilot/",code,".rda"))
    peptide_ids <- subset(get(paste("mapping",code,sep="_")),Protein==paste0(protein,"_HUMAN"))[["Isotope.Group.ID"]]
    peptide_exprs <- exprs(get(paste("peptide",code,sep="_")))
    peptide_exprs <- subset(peptide_exprs, rownames(peptide_exprs) %in% peptide_ids)
    peptide_dat <- dat %>%
                   left_join(data.frame(id=colnames(peptide_exprs),t(peptide_exprs))) %>%
                   filter(!is.na(FID) & grepl(code,id)) %>%
                   select(-id)
    peptides <- setdiff(names(peptide_dat),c("FID","IID","batch",ccovars))
    mod <- model.matrix(as.formula(paste0("~",paste(covars,collapse="+"))), data=peptide_dat)
    z <- sapply(names(peptide_dat[peptides]),
                function(col,verbose=FALSE)
                {
                  if (verbose) cat(names(peptide_dat[col]),col,"\n")
                  y <- gap::invnormal(peptide_dat[[col]])
                  l <- lm(y~mod[,-1])
                  r <- y-predict(l,na.action=na.pass)
                  scale(r)
                })
    colnames(z) <- names(peptide_dat[peptides])
    rownames(z) <- peptide_dat[["IID"]]
    peptide_lr <- data.frame(peptide_dat[c("FID","IID")],z)
    names(peptide_lr) <- gsub("^X([0-9])","\\1",names(peptide_lr))
    na_dat <- apply(peptide_dat[peptides],1,sum,na.rm=TRUE)
    if (verbose) write.table(filter(peptide_lr,na_dat!=0),
                             file=paste0(analysis,"/peptide/",protein,"/",code,".pheno"),
                             row.names=FALSE,quote=FALSE)
    peptide_lr
  }
  ZWK <- normalise_peptide(1)
  ZYQ <- normalise_peptide(2)
  UDP <- normalise_peptide(3)
  caprion <- bind_rows(ZWK,ZYQ,UDP)
  write.table(caprion,file=paste0(analysis,"/peptide/",protein,"/",protein,".pheno"),
              col.names=gsub("^X([0-9])","\\1",names(caprion)),row.names=FALSE,quote=FALSE)
  write.table(caprion,file=paste0(analysis,"/peptide/",protein,"/",protein,".mpheno"),
              col.names=FALSE,row.names=FALSE,quote=FALSE)
'

(
echo "#"'!'"/usr/bin/bash"
echo
echo "#SBATCH --account CARDIO-SL0-CPU"
echo "#SBATCH --partition cardio"
echo "#SBATCH --qos=cardio"
echo "#SBATCH --mem=28800"
echo "#SBATCH --time=12:00:00"
echo "#SBATCH --job-name=${protein}"
echo "#SBATCH --output=${analysis}/peptide/${protein}/${protein}.o"
echo "#SBATCH --error=${analysis}/peptide/${protein}/${protein}.e"
echo
echo "export protein=${protein}"
) > ${sbatch}
cat << 'EOL' >> ${sbatch}

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis

function fastLR()
{
  export batch=${1}
  export pheno=${analysis}/peptide/${protein}/${protein}.pheno
  export N=$(awk 'NR==1{print NF-2}' ${pheno})

  for col in $(seq ${N})
  do
  export peptide=$(awk 'NR==1{print $(col+2)}' col=${col} ${pheno})
  gcta-1.9 --mbgen ${pilot}/work/caprion.bgenlist \
           --sample ${pilot}/work/caprion.sample \
           --keep ${pilot}/work/caprion-${batch}.id \
           --fastGWA-lr \
           --pheno ${analysis}/peptide/${protein}/${protein}.mpheno --mpheno ${col} \
           --threads 10 \
           --out ${analysis}/peptide/${protein}/${protein}-${batch}-${peptide}

  gcta-1.9 --mbgen ${pilot}/work/caprion.bgenlist \
           --sample ${pilot}/work/caprion.sample \
           --keep ${pilot}/work/chrX-${batch}.id \
           --fastGWA-lr --model-only \
           --pheno ${analysis}/peptide/${protein}/${protein}.mpheno --mpheno ${col} \
           --threads 10 \
           --out ${analysis}/peptide/${protein}/${protein}-${batch}-${peptide}

  gcta-1.9 --bgen ${pilot}/work/chrX.bgen \
           --sample ${pilot}/work/chrX.sample \
           --keep ${pilot}/work/chrX-${batch}.id \
           --load-model ${analysis}/peptide/${protein}/${protein}-${batch}-${peptide}.fastGWA \
           --extract ${pilot}/work/chrX.snplist --geno 0.1 \
           --threads 10 \
           --out ${analysis}/peptide/${protein}/${protein}-${batch}-${peptide}-chrX
  gzip -f ${analysis}/peptide/${protein}/${protein}-${batch}-${peptide}*fastGWA
  done
}

fastLR 1
fastLR 2
fastLR 3
EOL

echo sbatch ${sbatch}
done
