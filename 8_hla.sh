#!/usr/bin/bash

export pilot=~/Caprion/pilot/data
export caprion=~/Caprion/analysis/work
export snp2hla=~/genetics/SNP2HLA_package_v1.0.3/SNP2HLA
export cookhla=${HPC_WORK}/CookHLA

function extract_hla()
{
  plink --bfile ${pilot}/merged_imputation --chr 6 --from-bp 25392021 --to-bp 33392022 \
        --make-bed --out ${caprion}/hla
}

function SNP2HLA()
{
  cd ${snp2hla}
# HapMap
  csh SNP2HLA.csh ${caprion}/hla ${snp2hla}/HM_CEU_REF ${caprion}/hla_IMPUTED plink 40000 100
# 1000Genomes
  csh SNP2HLA.csh ${caprion}/hla ${snp2hla}/1000G_REF.EUR.chr6.hg18.29mb-34mb.inT1DGC ${caprion}/hla_IMPUTED plink 40000 100
  plink --noweb --dosage ${caprion}/hla_IMPUTED.dosage noheader format=1 \
        --fam ${caprion}/hla_IMPUTED.fam \
        --linear --out ${caprion}/hla_IMPUTED.dosage.assoc
}

SNP2HLA

# CookHLA
# sbatch 8_hla.sb
# HIBAG
# R --no-save < 8_hla.R
