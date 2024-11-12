#!/usr/bin/bash

function sb()
# select from a list of peptide on duplicated proteins (1-706)
{
export analysis=~/Caprion/analysis
export sbatch=${analysis}/dup/dup-pgwas.sb
cat << 'EOL' > ${sbatch}
#!/usr/bin/bash

#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --job-name=_dup
#SBATCH --array=1-RUNS
#SBATCH --output=ANALYSIS/dup/slurm/_dup_%A_%a.o
#SBATCH --error=ANALYSIS/dup/slurm/_dup_%A_%a.e

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/htslib/1.20

export TMPDIR=${HPC_WORK}/work
export analysis=~/Caprion/analysis
export pheno=${analysis}/dup/ZWK.pheno

function fastLR()
{
  export batch=1
  export fastGWA=gcta-1.9
  export col=${SLURM_ARRAY_TASK_ID}
  export peptide=$(awk 'NR==1{print $(col+2)}' col=${col} ${pheno})
  export root=${analysis}/dup/ZWK
  ${fastGWA} --mbgen ${analysis}/bgen/caprion.bgenlist \
             --sample ${analysis}/bgen/caprion.sample \
             --extract ${analysis}/bgen/caprion.snplist \
             --keep ${analysis}/output/caprion-${batch}.id --geno 0.1 \
             --fastGWA-lr \
             --pheno ${root}.mpheno --mpheno ${col} \
             --threads 10 \
             --out ${root}-${batch}-${peptide}

  ${fastGWA} --mbgen ${analysis}/bgen/caprion.bgenlist \
             --sample ${analysis}/bgen/caprion.sample \
             --extract ${analysis}/bgen/caprion.snplist \
             --keep ${analysis}/output/chrX-${batch}.id --geno 0.1 \
             --fastGWA-lr --model-only \
             --pheno ${root}.mpheno --mpheno ${col} \
             --threads 10 \
             --out ${root}-${batch}-${peptide}-model

  ${fastGWA} --bgen ${analysis}/bgen/chrX.bgen \
             --sample ${analysis}/bgen/chrX.sample \
             --extract ${analysis}/bgen/chrX.snplist --geno 0.1 \
             --keep ${analysis}/output/chrX-${batch}.id \
             --load-model ${root}-${batch}-${peptide}-model.fastGWA \
             --threads 10 \
             --out ${root}-${batch}-${peptide}-chrX
  bgzip -f ${root}-${batch}-${peptide}.fastGWA
  tabix -f -S1 -s1 -b3 -e3 ${root}-${batch}-${peptide}.fastGWA.gz
  bgzip -f ${root}-${batch}-${peptide}-chrX.fastGWA
  tabix -f -S1 -s1 -b3 -e3 ${root}-${batch}-${peptide}-chrX.fastGWA.gz
}

fastLR 1
EOL

export analysis=~/Caprion/analysis
export pheno=${analysis}/dup/ZWK.pheno
export N=$(awk 'NR==1{print NF-2}' ${pheno})
sed -i "s|ANALYSIS|${analysis}|;s|RUNS|${N}|" ${sbatch}
sbatch ${sbatch}
}

function setup()
{
module load ceuadmin/R
Rscript -e '
    .libPaths()
    caprion <- "~/Caprion"
    analysis <- file.path(caprion,"analysis")
    load(file.path(analysis,"work","eSet.rda"))
    raw <- subset(raw_ZWK,grepl("\\|",raw_ZWK[[2]]))
    write.table(raw[1:6],file=file.path(analysis,"dup","dup.list"),col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
    suppressMessages(library(Biobase))
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyr))
    dup <- raw[-(2:6)]
    dupdat <- dup %>%
              tidyr::pivot_longer(cols = -Isotope.Group.ID,
                                  names_to = "ID",
                                  values_to = "Value") %>%
              dplyr::mutate(Value=if_else(Value==0,0,log2(Value))) %>%
              tidyr::pivot_wider(names_from = Isotope.Group.ID, values_from = Value)
    pcs <- paste0("ppc",1:3)
    ccovars <- c("sexPulse","agePulse",pcs,paste0("PC",1:20))
    dat <- read.delim(file.path(analysis,"output","caprion.pheno"))[c("FID","IID","caprion_id",ccovars)] %>%
           rename(id=caprion_id) %>%
           left_join(dupdat,by=c("id"="ID")) %>%
           filter(grepl("ZWK",id) & !is.na(FID))
    normalise_peptide <- function(batch,verbose=FALSE)
    {
      if (batch==1) covars <- setdiff(ccovars,pcs) else covars <- ccovars
      peptide_lr <- dat
      peptides <- setdiff(names(dat),c("FID","IID","id",ccovars))
      mod <- model.matrix(as.formula(paste0("~",paste(covars,collapse="+"))), data=dat[-c(1:3)])
      z <- sapply(names(dat[peptides]),
                  function(col)
                  {
                    y <- gap::invnormal(dat[[col]])
                    l <- lm(y~mod[,-1])
                    r <- y-predict(l,na.action=na.pass)
                    x <- gap::invnormal(r)
                    peptide_lr[col] <<- x
                    x
                  })
      names(peptide_lr[peptides]) <- gsub("^X([0-9])","\\1",names(peptide_lr[peptides]))
      peptide_lr[c("FID","IID",peptides)]
    }
    ZWK <- tryCatch(
      {
        normalise_peptide(1)
      },
      error = function(e) {
        cat("An error occurred:", conditionMessage(e), "\n")
        return(data.frame())
      },
      finally = {
        cat("Finished ZWK block.\n")
      }
    )
    write.table(ZWK,file=file.path(analysis,"dup","ZWK.pheno"),
                col.names=gsub("^X([0-9])","\\1",names(ZWK)),row.names=FALSE,quote=FALSE)
    write.table(ZWK,file=file.path(analysis,"dup","ZWK.mpheno"),
                col.names=FALSE,row.names=FALSE,quote=FALSE)
'
}

setup
sb
