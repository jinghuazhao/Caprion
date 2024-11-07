#!/usr/bin/bash

function sb()
# select from a list of duplicated proteins (either continuously in 1-706 or its subset)
{
  export index=${1}
  export protein=$(awk 'NR==ENVIRON["index"]{print $1}' ${analysis}/dup/dup.list)
  export sbatch=${analysis}/dup/dup-pgwas.sb

  Rscript -e '
    suppressMessages(library(Biobase))
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyr))
    .libPaths()
    caprion <- "~/Caprion"
    pilot <- file.path(caprion,"pilot")
    analysis <- file.path(caprion,"analysis")
    load(file.path(caprion,"pilot","ZWK.rda"))
    load(file.path(caprion,"analysis","work","eSet.rda"))
    dup <- subset(raw_ZWK[-(2:6)],grepl("\\|",raw_ZWK[[2]]))
    dupdat <- dup %>%
              tidyr::pivot_longer(cols = -Isotope.Group.ID,
                                  names_to = "ID",
                                  values_to = "Value") %>%
             tidyr::pivot_wider(names_from = Isotope.Group.ID, values_from = Value)
    pcs <- paste0("ppc",1:3)
    ccovars <- c("sexPulse","agePulse",pcs,paste0("PC",1:20))
    dat <- read.delim(file.path(analysis,"output","caprion.pheno"))[c("FID","IID","caprion_id",ccovars)] %>%
           rename(id=caprion_id) %>%
           filter(grepl("ZWK",id)) %>%
           left_join(dupdat,by=c("id"="ID"))
    isotope <- Sys.getenv("isotope")
    normalise_peptide <- function(batch,batches=c("ZWK","ZYQ","UDP"),verbose=FALSE)
    {
      if (batch==1) covars <- setdiff(ccovars,pcs) else covars <- ccovars
      peptide_ids <- dup[["Isotope.Group.ID"]]
      peptide_exprs <- dup[-1]
      peptide_dat <- dat %>%
                     select(-id)
      peptide_lr <- peptide_dat
      peptides <- setdiff(names(peptide_dat),c("FID","IID","batch",ccovars))
      mod <- model.matrix(as.formula(paste0("~",paste(covars,collapse="+"))), data=peptide_dat)
      z <- sapply(names(peptide_dat[peptides]),
                  function(col)
                  {
                    y <- gap::invnormal(peptide_dat[[col]])
                    l <- lm(y~mod[,-1])
                    r <- y-predict(l,na.action=na.pass)
                    x <- gap::invnormal(r)
                    peptide_lr[col] <<- x
                    x
                  })
      names(peptide_lr[peptides]) <- gsub("^X([0-9])","\\1",names(peptide_lr[peptides]))
      na_dat <- apply(peptide_dat[peptides],1,sum,na.rm=TRUE)
      peptide_table <- group_by(get(paste("mapping",code,sep="_")),Protein) %>%
                       summarise(N=n()) %>%
                       data.frame() %>%
                       filter(Protein!="-") %>%
                       select(N) %>%
                       table
      if (verbose) hist(peptide_table)
      if (verbose) write.table(filter(peptide_lr,na_dat!=0),
                               file=file.path(analysis,"dup","dup.pheno"),
                               row.names=FALSE,quote=FALSE)
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
    write.table(ZWK,file=file.path(analysis,"dup","dup.pheno"),
                col.names=gsub("^X([0-9])","\\1",names(caprion)),row.names=FALSE,quote=FALSE)
    write.table(ZWK,file=file.path(analysis,"dup","dup.mpheno"),
                col.names=FALSE,row.names=FALSE,quote=FALSE)
  '
  export pheno=${analysis}/dup/dup.pheno
  export N=$(awk 'NR==1{print NF-2}' ${pheno})
cat << 'EOL' > ${sbatch}
#!/usr/bin/bash

#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --job-name=PROTEIN
#SBATCH --array=1-RUNS
#SBATCH --output=ANALYSIS/peptide/PROTEIN/PROTEIN.o
#SBATCH --error=ANALYSIS/peptide/PROTEIN/PROTEIN.e

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load samtools/1.13/gcc/zwxn7ug3
module load ceuadmin/htslib/1.20

export isotope=ISOTOPE
export TMPDIR=${HPC_WORK}/work
export analysis=~/Caprion/analysis

function fastLR()
{
  export batch=${1}
  export fastGWA=gcta-1.9
  export col=${SLURM_ARRAY_TASK_ID}
  export peptide=$(awk 'NR==1{print $(col+2)}' col=${col} ${pheno})
  export root=${analysis}/peptide/${protein}/${protein}
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

sed -i "s|ANALYSIS|${analysis}|;s|ISOTOPE|${isotope}|g;s|RUNS|${N}|" ${sbatch}
sbatch ${sbatch}
}

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export varlist=${analysis}/dup/dup.list

module load ceuadmin/R

# all proteins:
while IFS=":" read -r isotope_index isotope; do
    export isotope_index
    export isotope
    echo ${isotope_index} ${isotope}
    export root=${analysis}/dup
  # sb ${isotope_index}
    export pheno=${root}/dup.pheno
    export N=$(awk "NR==1{print NF-2}" ${pheno})
done < <(cut -f1 ${varlist} | awk '{print NR":"$1}'
)
