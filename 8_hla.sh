#!/usr/bin/bash

export caprion=~/Caprion/analysis
export work=~/Caprion/pilot/work
export pilot=~/Caprion/pilot/data
export hatk=~/hpc-work/HATK/

function extract_hla()
# Region as used in the SCALLOP-INF project
{
  plink --bfile ${pilot}/merged_imputation --chr 6 --from-bp 25392021 --to-bp 33392022 \
        --keep ${work}/caprion.id2 \
        --make-bed --out ${caprion}/work/hla
}

extract_hla

# HLA imputation:
# SNP2HLA, CookHLA, HIBAG are now all part of the following SLURM script,
# sbatch 8_hla.sb
# HIBAG is specifically run through from the following R script,
# R --no-save < 8_hla.R
# CookHLA also uses SNP2HLA utility for making reference.
