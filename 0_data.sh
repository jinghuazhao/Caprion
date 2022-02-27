#!/usr/bin/bash

export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
export caprion=${HOME}/Caprion/pilot
export TMPDIR=${HPC_WORK}/work
export work=${caprion}/work
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
export phase=${2}
bgenix -g ${caprion}/${dir}/caprion.01.bgen -incl-rsids ${work}/variants.extract > ${work}/${phase}.bgen
bgenix -g ${work}/${phase}.bgen -index -clobber

sqlite3 ${work}/${phase}.bgen.bgi <<END
.tables
.separator "\t"
.header on
.output ${work}/metadata.txt
select * from metadata;
.output ${work}/Variant.txt
select * from Variant;
END

plink2 --bgen ${work}/${phase}.bgen ref-first --make-pfile --out ${work}/${phase}

set -eu

plink2 \
    --pfile ${work}/${phase} \
    --score ${PCA_loadings} \
            variance-standardize \
            cols=-scoreavgs,+scoresums \
            list-variants \
            header-read \
    --score-col-nums 3-12 \
    --read-freq ${PCA_af} \
    --out ${work}/${phase}

awk '
  {
    if (NR==1) $1="#FID IID"
    else $1=0" "$1
    print
  }' ${work}/${phase}.sscore > ${work}/${phase}-rsid.sscore
ln -sf ${work}/${phase}.sscore.vars ${work}/${phase}-rsid.sscore.vars

cat <(echo FID IID dummy | tr ' ' '\t') \
    <(sed '1d' ${work}/${phase}.psam | awk -v OFS='\t' '{print 0,$1,0}') > ${work}/${phase}-rsid.pheno

export ethnic_file=${HGI}/ethnic.txt
cat <(head -1  ${ethnic_file}) \
    <(cut -f1 ${work}/${phase}.psam | grep -f - ${ethnic_file}) > ${work}/${phase}-ethnic.txt
export covar_file=${HGI}/work/snpid.covars
cat <(head -1  ${covar_file}) \
    <(cut -f1 ${work}/${phase}.psam | grep -f - ${covar_file}) > ${work}/${phase}-covars.txt

Rscript ${PCA_projection}/plot_projected_pc.R \
        --sscore ${work}/${phase}-rsid.sscore \
        --phenotype-file ${work}/${phase}-rsid.pheno \
        --phenotype-col dummy \
        --covariate-file ${work}/${phase}-covars.txt \
        --pc-prefix PC \
        --pc-num 20 \
        --ancestry-file ${work}/${phase}-ethnic.txt \
        --ancestry-col ethnic \
        --study phase1 \
        --out ${work}/${phase}
}

pca_project data phase1
pca_project data2 phase2
pca_project data3 phase3

# qctool to combine bgen -- handling MAF<=0.01
# fine-mapping

# 196, 1488, 807 ==> 2491
# [data|data2|data3]/caprion.fam, caprion.bgen.bed/bim/fam, caprion.01.bgen/bgen.bgi
# [data|data2]/interval.sample

# ${caprion}/data/caprion.sample
# ${caprion}/data/phase2.sample
# ${caprion}/data3/protein.sample
