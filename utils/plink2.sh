#!/usr/bin/bash

export TMPDIR=/rds/user/$USER/hpc-work/

module load plink/2.00-alpha

function binary_ped()
{
seq 987 | \
parallel -C' ' '
  export col=$(cut -d" " -f {} caprion.uniprot); \
  for v in ${col} ${col}_invn;
  do
      echo ${v}
      plink2 \
             --bfile caprion.01 --fam caprion.fam \
             --glm hide-covar --input-missing-phenotype -9 --covar-variance-standardize \
             --pheno caprion.pheno --pheno-name ${v} --covar caprion.covar \
             --out work/${v}
      grep -v NA work/${v}.${v}.glm.linear | \
      gzip -f > plink2/${v}-plink2.gz
      rm work/${v}.${v}.glm.linear
  done
'
}

function bgen()
{
seq 987 | \
parallel -C' ' '
  export col=$(cut -d" " -f {} caprion.uniprot); \
  for v in ${col} ${col}_invn;
  do
      echo ${v}
      plink2 \
             --bgen caprion.01.bgen --sample caprion.sample \
             --glm hide-covar --input-missing-phenotype -9 --covar-variance-standardize \
             --pheno caprion.pheno --pheno-name ${v} --covar caprion.covar \
             --out work/${v}
      grep -v NA work/${v}.${v}.glm.linear | \
      gzip -f > bgen/${v}-plink2.gz
      rm work/${v}.${v}.glm.linear
  done
'
}

bgen
