#!/usr/bin/bash

export caprion=~/Caprion/pilot
export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37
export X=~/rds/projects/covid/ace2/interval_genetic_data/interval_imputed_data/
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
export TMPDIR=${HPC_WORK}/work

function X()
{
  awk '!/NA/{print $1"_"$1}' ${caprion}/work/caprion.id | \
  bcftools view -S - --force-samples ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -O z -o ${caprion}/work/INTERVAL-X.vcf.gz
  bcftools query -l ${caprion}/work/INTERVAL-X.vcf.gz | \
  tr '_' ' ' | \
  awk '{print $1}' | \
  bcftools reheader -s - ${caprion}/work/INTERVAL-X.vcf.gz -o ${caprion}/work/X.vcf.gz --threads 12
  bcftools index -tf ${caprion}/work/X.vcf.gz
  plink2 --vcf ${caprion}/work/X.vcf.gz --export bgen-1.2 bits=8 --double-id --dosage-erase-threshold 0.001 \
         --out ~/Caprion/analysis/work/X
  bgenix -g ~/Caprion/analysis/work/X.bgen -index -clobber
}
# bcftools annotate --set-id '%CHROM:%POS\_%REF\/%FIRST_ALT' ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -O z -o work/INTERVAL-X-vcf.gz

function fastGWAsetup()
{
# a sparse GRM from SNP data
  gcta-1.9 --bfile ${caprion}/work/caprion --make-grm --out ${caprion}/work/caprion
  gcta-1.9 --grm ${caprion}/work/caprion --make-bK-sparse 0.05 --out ${caprion}/work/caprion-spgrm
  echo ${caprion}/work/chr{1..22}.bgen | tr ' ' '\n' > ${caprion}/work/caprion.bgenlist
  (
    head -2 ${caprion}/work/chrX.sample
    sed '1,2d' ${caprion}/work/chrX.sample | \
    cut -d' ' -f1-3 | \
    join -a1 -e NA -o1.1,1.2,1.3,2.2 - ${caprion}/work/caprion.sex
  ) > ${caprion}/work/caprion.sample

  bcftools query -l ${caprion}/work/chrX.vcf.gz | awk '{print $1,$1}' > ${caprion}/work/chrX.idlist
  bcftools query -f "%ID\n" ${caprion}/work/chrX.vcf.gz > ${caprion}/work/chrX.snplist

# fastGWA mixed model
  sed -i '1d' ${caprion}/work/caprion.pheno
  sed -i '1d' ${caprion}/work/caprion-lr.pheno

  cut -f1,2 --complement ${caprion}/work/${caprion}.pheno | \
  head -1 | \
  tr '\t' '\n' > ${caprion}/work/caprion.varlist
}

# sbatch --export=ALL ${caprion}/5_pgwas.sb

function collect()
{
  ls ${caprion}/work/*.fastGWA.gz | \
  xargs -l basename | \
  sed 's/wes-//g;s/wgs-//g;s/.fastGWA.gz//g' | \
  grep -v chrX | \
  grep -v -f - ${caprion}/work/caprion.varlist > ${caprion}/work/caprion.lrlist

  parallel -j10 --env caprion -C' ' '
    export olink_protein={1}
    if [ -f ${caprion}/work/{2}-{1}.fastGWA.gz ] && [ -f ${caprion}/work/{2}-{1}-chrX.fastGWA.gz ]; then
       export out=${olink_protein}.txt.bgz
       (
         echo -e "CHR\tSNP\tPOS\tEFF_ALLELE\tOTHER_ALLELE\tN\tEFF_ALLELE_FREQ\tBETA\tSE\tP"
         gunzip -c ${caprion}/work/spa/{2}-{1}.fastGWA.gz | \
         sed "1d"
         gunzip -c ${caprion}/work/spa/{2}-{1}-chrX.fastGWA.gz | \
         sed "1d"
       ) | \
       awk -vOFS="\t" "{print \$2,\$1,\$3,\$6,\$4,\$5,\$7,\$8,\$9,\$10}" | \
       bgzip -f > ${caprion}/work/${out}
       tabix -f -S1 -s2 -b3 -e3 {caprion}/work/${out}
    fi
  ' ::: $(cat ${caprion}/work/caprion.varlist)
}

# Rscript -e 'library(HIBAG)'

function bgen()
{
  for chr in chr{1..22}
  do
     plink2 --vcf ${caprion}/work/${chr}.vcf.gz --export bgen-1.2 bits=8 --double-id --dosage-erase-threshold 0.001 \
            --out ${caprion}/work/${chr}
     bgenix -g ${caprion}/work/${chr}.bgen -index -clobber
  done
}

function tableMAF()
{
  seq 22 | \
  xargs -I {} cut -f15 ${ref}/impute_{}_interval.snpstats | \
  grep -v MAF | \
  sed '/^$/d' | \
  Rscript -e 'maf <- scan("stdin"); cutoffs <- cut(maf,c(0,0.001,0.01,0.05,0.1,0.2,0.3,0.4));table(cutoffs);table(cutoffs)/length(maf)'
}

# Read 87696888 items
# cutoffs
#    (0,0.001] (0.001,0.01]  (0.01,0.05]   (0.05,0.1]    (0.1,0.2]    (0.2,0.3]
#     63805560      9376909      2990789      1414448      1825305      1393653
#    (0.3,0.4]
#      1236541
# cutoffs
#    (0,0.001] (0.001,0.01]  (0.01,0.05]   (0.05,0.1]    (0.1,0.2]    (0.2,0.3]
#   0.72756926   0.10692408   0.03410371   0.01612883   0.02081379   0.01589170
#    (0.3,0.4]
#   0.01410017

# by chromosome
# seq 22 | xargs -I {} bash -c "cut -f15 ${ref}/impute_{}_interval.snpstats > {}.MAF"
