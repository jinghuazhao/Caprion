#!/usr/bin/bash

function sb()
# select from a list of all proteins (either continuously in 1-987 or its subset)
{
  export index=${1}
  export protein=$(awk 'NR==ENVIRON["index"]{print $1}' ${analysis}/output/caprion${suffix}.varlist)
  if [ ! -d ${analysis}/peptide/${protein} ]; then mkdir ${analysis}/peptide/${protein}; fi
  export sbatch=${analysis}/peptide/${protein}/${protein}-pgwas.sb

  Rscript -e '
    .libPaths()
    suppressMessages(library(Biobase))
    suppressMessages(library(dplyr))
    pilot <- Sys.getenv("pilot")
    analysis <- Sys.getenv("analysis")
    protein <- Sys.getenv("protein")
    pcs <- paste0("ppc",1:3)
    ccovars <- c("sexPulse","agePulse",pcs,paste0("PC",1:20))
    dat <- read.delim(file.path(analysis,"output","caprion.pheno"))[c("FID","IID","caprion_id",ccovars)] %>%
           rename(id=caprion_id)
    normalise_peptide <- function(batch,batches=c("ZWK","ZYQ","UDP"),verbose=FALSE)
    {
      if (batch==1) covars <- setdiff(ccovars,pcs) else covars <- ccovars
      code <- batches[batch]
      load(paste0(pilot,"/",code,".rda"))
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
                    gap::invnormal(r)
                  })
      colnames(z) <- names(peptide_dat[peptides])
      rownames(z) <- peptide_dat[["IID"]]
      peptide_lr <- data.frame(peptide_dat[c("FID","IID")],z)
      names(peptide_lr) <- gsub("^X([0-9])","\\1",names(peptide_lr))
      na_dat <- apply(peptide_dat[peptides],1,sum,na.rm=TRUE)
      peptide_table <- group_by(get(paste("mapping",code,sep="_")),Protein) %>%
                       summarise(N=n()) %>%
                       data.frame() %>%
                       filter(Protein!="-") %>%
                       select(N) %>%
                       table
      if (verbose) hist(peptide_table)
      if (verbose) write.table(filter(peptide_lr,na_dat!=0),
                               file=paste0(analysis,"/peptide/",protein,"/",code,".pheno"),
                               row.names=FALSE,quote=FALSE)
      peptide_lr
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
    ZYQ <- tryCatch(
      {
        normalise_peptide(2)
      },
      error = function(e) {
        cat("An error occurred:", conditionMessage(e), "\n")
        return(data.frame())
      },
      finally = {
        cat("Finished ZYQ block.\n")
      }
    )
    UDP <- tryCatch(
      {
        normalise_peptide(3)
      },
      error = function(e) {
        cat("An error occurred:", conditionMessage(e), "\n")
        return(data.frame())
      },
      finally = {
        cat("Finished UDP block.\n")
      }
    )
    caprion <- bind_rows(ZWK,ZYQ,UDP)
    write.table(caprion,file=paste0(analysis,"/peptide/",protein,"/",protein,".pheno"),
                col.names=gsub("^X([0-9])","\\1",names(caprion)),row.names=FALSE,quote=FALSE)
    write.table(caprion,file=paste0(analysis,"/peptide/",protein,"/",protein,".mpheno"),
                col.names=FALSE,row.names=FALSE,quote=FALSE)
  '
  export pheno=${analysis}/peptide/${protein}/${protein}.pheno
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

export protein=PROTEIN
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
  bgzip -f ${root}-${batch}-${peptide}-chrX.fastGWA
}

fastLR 1
fastLR 2
fastLR 3
EOL

sed -i "s|ANALYSIS|${analysis}|;s|PROTEIN|${protein}|g;s|RUNS|${N}|" ${sbatch}
sbatch ${sbatch}
}

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export signals=${analysis}/work/caprion${suffix}.signals
export varlist=${analysis}/output/caprion${suffix}.varlist

if [ "$(uname -n | sed 's/-[0-9]*$//')" == "login-q" ]; then
   module load ceuadmin/R/4.4.0-icelake
else
   module load ceuadmin/R
fi

# all proteins:
xargs -n 2 < ${analysis}/peptide_progs/benchmark2.names | \
grep -n -f ${analysis}/peptide_progs/benchmark2.names -w ${varlist} | \
while IFS=":" read -r protein_index protein; do
    export protein_index
    export protein
    echo ${protein_index} ${protein}
    export root=${analysis}/peptide/${protein}
    sb ${protein_index}
    export pheno=${root}/${protein}.pheno
    export N=$(awk "NR==1{print NF-2}" ${pheno})
done
