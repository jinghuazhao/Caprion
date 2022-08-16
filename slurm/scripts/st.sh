#!/usr/bin/bash

module load miniconda3/4.5.1
export csd3path=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop
source activate ${csd3path}/miniconda37
snakemake -c --profile slurm -s slurm/rules/cojo.yaml --unlock
scancel -u jhz22
snakemake -c --profile slurm -s slurm/rules/cojo.yaml

# snakemake -c --cluster-config slurm/config.yaml -s slurm/rules/cojo.yaml -d slurm # temptating syntax

