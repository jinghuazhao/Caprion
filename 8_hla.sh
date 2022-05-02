#!/usr/bin/bash

export pilot=~/Caprion/pilot/data
export caprion=~/Caprion/analysis/work

function extract_hla()
# Region as used in the SCALLOP-INF project
{
  plink --bfile ${pilot}/merged_imputation --chr 6 --from-bp 25392021 --to-bp 33392022 \
        --make-bed --out ${caprion}/hla
}

# HLA imputation:
# SNP2HLA, CookHLA, HIBAG are now all part of the following SLURM script,
# sbatch 8_hla.sb
# HIBAG is specifically run through from the following R script,
# R --no-save < 8_hla.R
# CookHLA also uses SNP2HLA utility for making reference.
