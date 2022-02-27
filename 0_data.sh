#!/usr/bin/bash

export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
export caprion=${HOME}/Caprion/pilot
export TMPDIR=${HPC_WORK}/work
export dir=~/Caprion/pilot/work
export HGI=~/COVID-19/HGI
export PCA_projection=${HGI}/pca_projection
export PCA_loadings=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.tsv
export PCA_af=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.afreq
export sscore=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz

module load plink/2.00-alpha ceuadmin/stata

cut -f1 ${PCA_loadings} | tail -n +2 > ${dir}/variants.extract
bgenix -g ${caprion}/data/caprion.01.bgen -incl-rsids ${dir}/variants.extract > ${dir}/caprion.bgen
bgenix -g ${dir}/caprion.bgen -index -clobber

sqlite3 ${dir}/caprion.bgen.bgi <<END
.tables
.separator "\t"
.header on
.output ${dir}/metadata.txt
select * from metadata;
.output ${dir}/Variant.txt
select * from Variant;
END

plink2 --bgen ${dir}/caprion.bgen ref-first --make-pfile --out ${dir}/caprion

set -eu

plink2 \
    --pfile ${dir}/caprion \
    --score ${PCA_loadings} \
            variance-standardize \
            cols=-scoreavgs,+scoresums \
            list-variants \
            header-read \
    --score-col-nums 3-12 \
    --read-freq ${PCA_af} \
    --out ${dir}/caprion

awk '
  {
    if (NR==1) $1="#FID IID"
    else $1=0" "$1
    print
  }' ${dir}/caprion.sscore > ${dir}/snpid.sscore
ln -sf ${dir}/caprion.sscore.vars ${dir}/snpid.sscore.vars

stata -b do ${HGI}/ethnic.do

Rscript ${PCA_projection}/plot_projected_pc.R \
        --sscore ${dir}/snpid.sscore \
        --phenotype-file ${dir}/snpid.pheno \
        --phenotype-col SARS_CoV \
        --covariate-file ${dir}/snpid.covars \
        --pc-prefix PC \
        --pc-num 20 \
        --ancestry-file ethnic.txt \
        --ancestry-col ethnic \
        --study INTERVAL \
        --out ${dir}/INTERVAL

project_pc

# qctool to combine bgen, extract ids
# PCA projection
# fine-mapping

# 196, 1488, 807 ==> 2491
# [data|data2|data3]/caprion.fam, caprion.bgen.bed/bim/fam, caprion.01.bgen/bgen.bgi
# [data|data2]/interval.sample

# ${caprion}/data/caprion.sample
# ${caprion}/data/phase2.sample
# ${caprion}/data3/protein.sample
