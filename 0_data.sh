#!/usr/bin/bash

export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
export caprion=${HOME}/Caprion/pilot
export TMPDIR=${HPC_WORK}/work
export dir=~/Caprion/pilot/work
export HGI=~/COVID-19/HGI
export PCA_projection=${HGI}/pca_projection
export PCA_loadings=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.tsv
export PCA_af=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.afreq
export sscore=${PCA_projection}/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz

module load plink/2.00-alpha ceuadmin/stata

function extract_data()
{
  cut -f1 ${PCA_loadings} | tail -n +2 > ${dir}/variants.extract
  (
    cat ${dir}/variants.extract
    awk '{split($1,a,":");print "chr"a[1]":"a[2]"_"a[4]"_"a[3]}' ${dir}/variants.extract
  ) > ${dir}/variants.extract2

  cp ${caprion}/data/caprion.01.bgen.bgi ${TMPDIR}
  seq 22 | \
  parallel -C' ' --env prefix --env dir '
    ln -sf ${caprion}/data/caprion.01.bgen ${TMPDIR}/caprion-{}.bgen
    python ${HGI}/update_bgi.py --bgi ${TMPDIR}/caprion-{}.bgen.bgi
    bgenix -g ${TMPDIR}/INTERVAL-{}.bgen -incl-rsids variants.extract2 > ${dir}/caprion-{}.bgen
  '
}

function twist()
{
  cat-bgen -g $(echo ${dir}/caprion-{1..22}.bgen) -og caprion.bgen -clobber
  bgenix -g caprion.bgen -index -clobber
  export csvfile=INTERVAL.csv
  python update_bgi.py --bgi caprion.bgen.bgi
  plink2 --bgen caprion.bgen ref-first --make-pfile --out caprion
  cp caprion.p??? ${dir}
  (
    head -1 caprion.pvar
    paste <(sed '1d' caprion.pvar | cut -f1,2) \
          <(sed '1d' INTERVAL.csv | cut -d, -f3) | \
    paste - <(sed '1d' caprion.pvar | cut -f4,5)
  ) > ${dir}/caprion.pvar
}

function project_pc()
{
  set -eu
  PCA_LOADINGS="${PCA_loadings}"
  PCA_AF="${PCA_af}"
  PFILE="${prefix}/work/INTERVAL"
  BFILE=""

  plink2 \
    --pfile INTERVAL \
    --score ${PCA_LOADINGS} \
            variance-standardize \
            cols=-scoreavgs,+scoresums \
            list-variants \
            header-read \
    --score-col-nums 3-12 \
    --read-freq ${PCA_af} \
    --out ${dir}/INTERVAL

  awk '
  {
    if (NR==1) $1="#FID IID"
    else $1=0" "$1
    print
  }' ${INTERVAL}.sscore > ${dir}/snpid.sscore
  ln -sf ${dir}/INTERVAL.sscore.vars ${dir}/snpid.sscore.vars

  stata -b do ${HGI}/ethnic.do

  Rscript ${PCA_projection}/plot_projected_pc.R \
        --sscore ${prefix}/work/snpid.sscore \
        --phenotype-file ${dir}/snpid.pheno \
        --phenotype-col SARS_CoV \
        --covariate-file ${dir}snpid.covars \
        --pc-prefix PC \
        --pc-num 20 \
        --ancestry-file ethnic.txt \
        --ancestry-col ethnic \
        --study INTERVAL \
        --out ${dir}/INTERVAL
}

extract_data;twist;project_pc

# qctool to combine bgen, extract ids
# PCA projection
# fine-mapping

# 196, 1488, 807 ==> 2491
# [data|data2|data3]/caprion.fam, caprion.bgen.bed/bim/fam, caprion.01.bgen/bgen.bgi
# [data|data2]/interval.sample

# ${caprion}/data/caprion.sample
# ${caprion}/data/phase2.sample
# ${caprion}/data3/protein.sample
