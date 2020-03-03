#!/usr/bin/bash

export TMPDIR=/rds/user/$USER/hpc-work/

module load plink/2.00-alpha

seq 987 | \
parallel -C' ' '
  export col=$(cut -d" " -f {} caprion.uniprot); \
  for v in ${col} ${col}_inv;
  do
      echo ${v}
      plink2 \
             --bfile caprion.01 --fam caprion.fam \
             --glm hide-covar --input-missing-phenotype -9 --covar-variance-standardize \
             --pheno caprion.pheno --pheno-name ${col} --covar caprion.covar \
             --out work/${v}
      grep -v NA work/${v}.${col}.glm.linear | \
      gzip -f > plink2/${v}-plink2.gz
      rm work/${v}.${col}.glm.linear
  done
'
