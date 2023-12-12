#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr

function sb()
{
  export index=${1}
# all proteins
  export protein=$(awk 'NR==ENVIRON["index"]{print $1}' ${pilot}/work/caprion.varlist)
# only those with pQTLs
# export protein=$(awk 'NR>1{print $1}' ${analysis}/work/caprion${suffix}.signals | sort -k1,1 | uniq | awk 'NR==ENVIRON["index"]')
  if [ ! -d ${analysis}/peptide/${protein} ]; then mkdir ${analysis}/peptide/${protein}; fi
  export sbatch=${analysis}/peptide/${protein}/${protein}-pgwas.sb

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
    ZWK <- normalise_peptide(1)
    ZYQ <- normalise_peptide(2)
    UDP <- normalise_peptide(3)
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
#SBATCH --partition cclake
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --job-name=PROTEIN
#SBATCH --array=1-RUNS
#SBATCH --output=ANALYSIS/peptide/PROTEIN/PROTEIN.o
#SBATCH --error=ANALYSIS/peptide/PROTEIN/PROTEIN.e

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-ccl

export protein=PROTEIN
export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis

function fastLR()
{
  export batch=${1}
  export fastGWA=gcta-1.9

  export col=${SLURM_ARRAY_TASK_ID}
  export peptide=$(awk 'NR==1{print $(col+2)}' col=${col} ${pheno})
  export root=${analysis}/peptide/${protein}/${protein}
  ${fastGWA} --mbgen ${pilot}/work/caprion.bgenlist \
             --sample ${pilot}/work/caprion.sample \
             --keep ${pilot}/work/caprion-${batch}.id --geno 0.1 --maf 0.001 \
             --fastGWA-lr \
             --pheno ${root}.mpheno --mpheno ${col} \
             --threads 10 \
             --out ${root}-${batch}-${peptide}

  ${fastGWA} --mbgen ${pilot}/work/caprion.bgenlist \
             --sample ${pilot}/work/caprion.sample \
             --keep ${pilot}/work/chrX-${batch}.id --geno 0.1 --maf 0.001 \
             --fastGWA-lr --model-only \
             --pheno ${root}.mpheno --mpheno ${col} \
             --threads 10 \
             --out ${root}-${batch}-${peptide}-model

  ${fastGWA} --bgen ${pilot}/work/chrX.bgen \
             --sample ${pilot}/work/chrX.sample \
             --keep ${pilot}/work/chrX-${batch}.id \
             --load-model ${root}-${batch}-${peptide}-model.fastGWA \
             --extract ${pilot}/work/chrX.snplist --geno 0.1 \
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

# for i in $(seq 987); do sb ${i}; done

for i in $(grep error ${analysis}/peptide/*/*.e | sed 's|/|\t|g' | cut -f7 | grep -n -w -f - ${pilot}/work/caprion.varlist | cut -d':' -f1)
do
  sb ${i}
done
