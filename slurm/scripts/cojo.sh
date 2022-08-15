#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export bfile=${1}
export p=${2}
export r=${3}
export pr=${p}-${r}
export chr=${4}
export pos=${5}
export flanking=1e6
export start=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{start=pos-flanking;if(start<0) start=0;print start}')
export end=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{print pos+flanking}')

function cojo()
{
  cat <(echo SNP A1 A2 freq b se p N) \
      <(gunzip -c METAL/sentinels/${p}.p.gz | awk '{print $3,toupper($4),toupper($5),$6,$10,$11,10^$12,$18}') > work/${p}.ma
  cut -d' ' -f1 work/${p}.ma | sed '1d' > work/${p}.rsid
  # single significant variant
  if [ $(wc -l work/${p}.rsid) -eq 1 ]; then return 0; fi
  plink2 --bgen data/chr${chr}.bgen ref-last --sample data/caprion.sample \
         --extract work/${p}.rsid --indep-pairwise 1000kb 1 0.1 \
         --out results/${pr}.prune
  # only one variant left
  if [ $(wc -l work/${pr}.prune.prune.in) -eq 1 ]; then return 0; fi
  if [ $(grep -w ${r} results/${pr}.prune.prune.in | wc -l) -eq 0 ]; then
     export i=$(grep -w -f results/${pr}.prune.prune.in ${bfile}.bim | \
                awk -vpos=${pos} 'function abs(x) {if (x<0) return -x; else return x;} {print $1, $2, $4, abs($4-pos)}' | \
                sort -r -k4,4n | \
                awk 'NR==1 {print $2}' \
               )
     sed -i 's/'"$i"'/'"$r"'/g' results/${pr}.prune.prune.in
  fi
  (
    if [ ${chr} -eq 19 ]; then
       sort results/${pr}.prune.prune.in | join -v1 - ${INF}/work/NLRP2
    else
       sort results/${pr}.prune.prune.in
    fi
  ) > results/${pr}.prune
# rm results/${pr}.prune.prune.in results/${pr}.prune.prune.out work/${p}.rsid
  plink2 --bgen data/chr${chr}.bgen ref-last \
         --sample data/caprion.sample \
         --extract results/${pr}.prune \
         --recode A include-alt \
         --out results/${pr}.dosage
  plink2 --bgen data/chr${chr}.bgen ref-last --sample data/caprion.sample \
         --extract results/${pr}.prune --export ind-major-bed --out work/${pr}
  gcta-1.9 --bfile work/${pr} \
           --cojo-file work/${p}.ma --extract-region-snp ${r} 1000 --maf 0.01 \
           --cojo-slct \
           --cojo-p 5e-8 \
           --cojo-collinear 0.9 \
           --out results/${pr}.gcta
  Rscript -e '
    p <- Sys.getenv("p")
    pr <- Sys.getenv("pr")
    r <- scan(paste0("results/",pr,".prune"),"")
    suppressMessages(library(dplyr))
    raw <- read.delim(paste0("results/",pr,".dosage.raw")) %>% select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)
    load("reports/edata_batch.rda")
    load("data/pheno.rda")
    edata <- with(edata_batch,edata)
    edata <- left_join(subset(id,!is.na(IID)),data.frame(caprion_id=rownames(edata),edata)) %>%
             select(IID,p) %>%
             left_join(raw) %>%
             select(-IID)
    intercept_only <- lm(edata[[p]] ~ 1, data=edata)
    all <- lm(edata[[p]] ~ ., data=edata)
    both <- step(intercept_only, direction='both', scope=formula(all), trace=0)
    both$anova
    summary(both)
    chunks <- split(r, ceiling(seq_along(r)/snakemake@params[["chunksize"]]))
  '
}

cojo

function setup()
{
  mkdir data
  cd data
  cat ~/Caprion/pilot/work/caprion.bgenlist | xargs -l -I {} ln -s {}
  ln -s ~/Caprion/pilot/work/caprion.sample
  save(pheno,id,file="data/pheno.rda")
  ln -s /home/jhz22/Caprion/pilot/work/caprion.pheno
  Rscript -e '
  pheno <- read.delim("caprion.pheno")
  id <- pheno[c("IID","caprion_id")]
  save(pheno,id,file="pheno.rda")
  '
  cd -
# https://www.statology.org/stepwise-regression-r/
}
