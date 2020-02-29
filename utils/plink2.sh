#!/usr/bin/bash

export TMPDIR=/rds/user/$USER/hpc-work/

module load plink/2.00-alpha

plink --bfile caprion.bgen --maf 0.01 --make-bed --out caprion.01
awk '{$1=$2};1' caprion.01.fam > caprion.fam
seq 609 | \
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
