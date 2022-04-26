#!/usr/bin/bash

export caprion=~/Caprion/pilot
export interval=${HPC_WORK}/data/interval
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
  bcftools query -l ${caprion}/work/X.vcf.gz | awk '{print $1,$1}' > ${caprion}/work/chrX.idlist
  bcftools query -f "%ID\n" ${caprion}/work/X.vcf.gz > ${caprion}/work/chrX.snplist
  plink2 --vcf ${caprion}/work/X.vcf.gz --export bgen-1.2 bits=8 --double-id --dosage-erase-threshold 0.001 \
         --set-missing-var-ids @:#_\$r_\$a --new-id-max-allele-len 680 \
         --out ${caprion}/work/chrX
  bgenix -g ${caprion}/work/chrX.bgen -index -clobber
  cut -d' ' -f1 ${caprion}/work/caprion-1.id | grep -f - ${caprion}/work/chrX.idlist > ${caprion}/work/chrX-1.id
  cut -d' ' -f1 ${caprion}/work/caprion-2.id | grep -f - ${caprion}/work/chrX.idlist > ${caprion}/work/chrX-2.id
  cut -d' ' -f1 ${caprion}/work/caprion-3.id | grep -f - ${caprion}/work/chrX.idlist > ${caprion}/work/chrX-3.id
}

function fastGWAsetup()
{
  paste -d' ' ${caprion}/work/caprion.id ${caprion}/work/caprion.id | grep -v NA > ${caprion}/work/caprion.id2
# echo ${caprion}/work/chr{1..22}.bgen | tr ' ' '\n' > ${caprion}/work/caprion.bgenlist
  seq 22 | \
  xargs -I {} echo ${caprion}/work/chr{}.bgen > ${caprion}/work/caprion.bgenlist
  cat <(head -2 ${caprion}/data/merged_imputation.sample | awk '{if(NR==1) {print $0, "sex"} else {print $0, "D"}}') \
      <(grep -f ${caprion}/work/caprion.id -w ${caprion}/data/merged_imputation.sample | \
        cut -d' ' -f1-2 | join - ${caprion}/data/merged_imputation.missing | \
        awk '{print $0, "NA"}') > ${caprion}/work/caprion.sample
  gcta-1.9 --bfile ${caprion}/data/merged_imputation --keep ${caprion}/work/caprion.id2 --make-grm --out ${caprion}/work/caprion --threads 10
# a sparse GRM from SNP data
  gcta-1.9 --grm ${caprion}/work/caprion --make-bK-sparse 0.05 --out ${caprion}/work/caprion-spgrm --threads 10
# fastGWA mixed model
  cut -f1,2 --complement ${caprion}/work/caprion.pheno | \
  head -1 | \
  tr '\t' '\n' > ${caprion}/work/caprion.varlist
  sed -i '1d' ${caprion}/work/caprion-1.pheno
  sed -i '1d' ${caprion}/work/caprion-2.pheno
  sed -i '1d' ${caprion}/work/caprion-3.pheno
}
# slow but unnecessary
# gcta-1.9 --mbgen ${caprion}/work/caprion.bgenlist --sample ${caprion}/work/caprion.sample --make-grm --out ${caprion}/work/caprion --threads 10
# sed -i '1d' ${caprion}/work/caprion.pheno

# sbatch --export=ALL ${caprion}/5_pgwas.sb

if [ ! -d ~/Caprion/analysis/METAL/qqmanhattan ]; then mkdir ~/Caprion/analysis/METAL/qqmanhattan; fi

function lrlist()
{
  ls ${caprion}/work/*.fastGWA.gz | \
  xargs -l basename | \
  sed 's/.fastGWA.gz//g' | \
  grep -v chrX | \
  grep -v -f - ${caprion}/work/caprion.varlist > ${caprion}/work/caprion.lrlist
}

function collect()
{
  parallel -j10 --env caprion -C' ' '
    export caprion_protein={1}
    if [ -f ${caprion}/work/{2}-{1}.fastGWA.gz ] && [ -f ${caprion}/work/{2}-{1}-chrX.fastGWA.gz ]; then
       export out=${caprion_protein}.txt.bgz
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

function tableMAF()
{
  seq 22 | \
  xargs -I {} cut -f15 ${ref}/impute_{}_interval.snpstats | \
  grep -v MAF | \
  sed '/^$/d' | \
  Rscript -e '
    maf <- scan("stdin");
    options(width=200)
    common_cutoffs <- cut(maf,c(0,0.01,0.05,0.5));
    table(common_cutoffs);
    table(common_cutoffs)/length(maf);
    full_cutoffs <- cut(maf,c(0,0.0001,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5));
    table(full_cutoffs);
    table(full_cutoffs)/length(maf)
  '
}

# Read 87696888 items
# common_cutoffs
#    (0,0.01] (0.01,0.05]  (0.05,0.5]
#    73182469     2990789     7038912
# common_cutoffs
#    (0,0.01] (0.01,0.05]  (0.05,0.5]
#  0.83449334  0.03410371  0.08026410
# full_cutoffs
#     (0,0.0001] (0.0001,0.001]   (0.001,0.01]    (0.01,0.05]     (0.05,0.1]      (0.1,0.2]      (0.2,0.3]      (0.3,0.4]      (0.4,0.5]
#       38883052       24922508        9376909        2990789        1414448        1825305        1393653        1236541        1168965
# full_cutoffs
#     (0,0.0001] (0.0001,0.001]   (0.001,0.01]    (0.01,0.05]     (0.05,0.1]      (0.1,0.2]      (0.2,0.3]      (0.3,0.4]      (0.4,0.5]
#     0.44338007     0.28418919     0.10692408     0.03410371     0.01612883     0.02081379     0.01589170     0.01410017     0.01332961

# by chromosome
# seq 22 | xargs -I {} bash -c "cut -f15 ${ref}/impute_{}_interval.snpstats > {}.MAF"