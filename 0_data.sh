#!/usr/bin/bash

export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
export caprion=${HOME}/Caprion
export pilot=${caprion}/pilot
export TMPDIR=${HPC_WORK}/work
export work=${pilot}/work
export HGI=~/COVID-19/HGI
export PCA_projection=${HGI}/pca_projection
export PCA_loadings=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.tsv
export PCA_af=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.afreq
export sscore=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz

module load plink/2.00-alpha ceuadmin/stata

cut -f1 ${PCA_loadings} | tail -n +2 > ${work}/variants.extract

function pca_project()
# PCA projection
{
export dir=${1}
export batch=${2}
bgenix -g ${pilot}/${dir}/caprion.01.bgen -incl-rsids ${work}/variants.extract > ${work}/${batch}.bgen
bgenix -g ${work}/${batch}.bgen -index -clobber

sqlite3 ${work}/${batch}.bgen.bgi <<END
.tables
.separator "\t"
.header on
.output ${work}/metadata.txt
select * from metadata;
.output ${work}/Variant.txt
select * from Variant;
END

plink2 --bgen ${work}/${batch}.bgen ref-first --make-pfile --out ${work}/${batch}

set -eu

plink2 \
    --pfile ${work}/${batch} \
    --score ${PCA_loadings} \
            variance-standardize \
            cols=-scoreavgs,+scoresums \
            list-variants \
            header-read \
    --score-col-nums 3-12 \
    --read-freq ${PCA_af} \
    --out ${work}/${batch}

awk '
  {
    if (NR==1) $1="#FID IID"
    else $1=0" "$1
    print
  }' ${work}/${batch}.sscore > ${work}/${batch}-rsid.sscore
ln -sf ${work}/${batch}.sscore.vars ${work}/${batch}-rsid.sscore.vars

cat <(echo FID IID dummy | tr ' ' '\t') \
    <(sed '1d' ${work}/${batch}.psam | awk -v OFS='\t' '{print 0,$1,0}') > ${work}/${batch}-rsid.pheno

export ethnic_file=${HGI}/ethnic.txt
cat <(head -1  ${ethnic_file}) \
    <(cut -f1 ${work}/${batch}.psam | grep -f - ${ethnic_file}) > ${work}/${batch}-ethnic.txt
export covar_file=${HGI}/work/snpid.covars
cat <(head -1  ${covar_file}) \
    <(cut -f1 ${work}/${batch}.psam | grep -f - ${covar_file}) > ${work}/${batch}-covars.txt

cut -f1 ${work}/${batch}.psam | \
grep -f - ${caprion}/INTERVAL_OmicsMap_20210915.csv | \
cut -d, -f1 | \
grep -f - ${caprion}/INTERVALdata_15SEP2021.csv | \
cut -d, -f1,7

Rscript ${PCA_projection}/plot_projected_pc.R \
        --sscore ${work}/${batch}-rsid.sscore \
        --phenotype-file ${work}/${batch}-rsid.pheno \
        --phenotype-col dummy \
        --covariate-file ${work}/${batch}-covars.txt \
        --pc-prefix PC \
        --pc-num 20 \
        --ancestry-file ${work}/${batch}-ethnic.txt \
        --ancestry-col ethnic \
        --study ${batch} \
        --out ${work}/${batch}
}

pca_project data batch1
pca_project data2 batch2
pca_project data3 batch3

# data handling with MAF<=0.01
# [data|data2|data3]/caprion.fam, caprion.bgen.bed/bim/fam, caprion.01.bgen/bgen.bgi
# [data|data2]/interval.sample, [data3]/[dr|peptide|protein].sample

# expr 196 + 1488 + 807
# ==> 2491

# ${pilot}/data/caprion.sample
# ${pilot}/data/phase2.sample
# ${pilot}/data3/protein.sample
# qctool to combine bgen fine-mapping?
