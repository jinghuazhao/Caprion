#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export bfile=${1}
export p=${2}
export r=${3}
export pr=${p}-${r}
export chr=${4}
export pos=${5}

function cojo()
{
  cat <(echo SNP A1 A2 freq b se p N) \
      <(gunzip -c METAL3/sentinels/${p}.p.gz | awk '{print $3,toupper($4),toupper($5),$6,$10,$11,10^$12,$18}') > work/${pr}.ma
  cut -d' ' -f1 work/${pr}.ma | sed '1d' > work/${pr}.rsid
  # singleton variant
  if [ $(wc -l work/${pr}.rsid | cut -d' ' -f1) -eq 1 ]; then return 0; fi
  plink2 --bgen data/chr${chr}.bgen ref-unknown --sample data/caprion.sample \
         --extract work/${pr}.rsid --export ind-major-bed --out work/${pr}
  gcta-1.9 --bfile work/${pr} \
           --cojo-file work/${pr}.ma --maf 0.01 --diff-freq 1 \
           --cojo-slct \
           --cojo-p 5e-8 \
           --cojo-collinear 0.9 \
           --out results/${pr}.gcta
  # for flanking window adding --extract-region-snp ${r} 1000
  if [ ! -f results/${pr}.gcta.jma.cojo ]; then return 0; fi
  sed '1d' results/${pr}.gcta.jma.cojo | cut -f2 > results/${pr}.jma
  plink --bfile work/${pr} \
        --extract results/${pr}.jma \
        --r square \
        --out results/${pr}
  plink2 --bgen data/chr${chr}.bgen ref-unknown \
         --sample data/caprion.sample \
         --extract results/${pr}.jma \
         --recode A include-alt \
         --out results/${pr}.dosage
  R --no-save < slurm/scripts/cojo.R > results/${pr}.lm.log
  rm work/${pr}.ma work/${pr}.rsid results/${pr}.jma
}

cojo
