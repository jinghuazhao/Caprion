#!/usr/bin/bash

# genotype
export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/
# location of PCs
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/

# Caprion working directories
export caprion=${HOME}/Caprion/pilot
export TMPDIR=${HPC_WORK}/work
export dir=~/Caprion/pilot/work

export PCA_projection=~/COVID-19/HGI/pca_projection
export PCA_loadings=hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.tsv
export PCA_af=hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.afreq
export sscore=hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz

module load plink/2.00-alpha ceuadmin/stata

function extract_data()
{
  cut -f1 ${PCA_projection}/${PCA_loadings} | tail -n +2 > ${dir}/variants.extract
  (
    cat ${dir}/variants.extract
    awk '{split($1,a,":");print "chr"a[1]":"a[2]"_"a[4]"_"a[3]}' ${dir}/variants.extract
  ) > ${dir}/variants.extract2

  cp ${caprion}/data/caprion.01.bgen.bgi ${TMPDIR}
  seq 22 | \
  parallel -C' ' --env prefix --env dir '
    ln -sf ${caprion}/data/caprion.01.bgen ${TMPDIR}/caprion-{}.bgen
    python update_bgi.py --bgi ${TMPDIR}/caprion-{}.bgen.bgi
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
# Metadata
  STUDY_NAME="INTERVAL"
  ANALYST_LAST_NAME="ZHAO"
  DATE="$(date +'%Y%m%d')"
  OUTNAME="$dir}/${STUDY_NAME}.${ANALYST_LAST_NAME}.${DATE}"
# Location of downloaded input files
  PCA_LOADINGS="${PCA_projection}/${PCA_loadings}"
  PCA_AF="${PCA_projection}/${PCA_af}"
# Location of imputed genotype files
# [Recommended]
# PLINK 2 binary format: a prefix (with directories) of .pgen/.pvar/.psam files
  PFILE="${prefix}/work/INTERVAL"
# [Acceptable]
# PLINK 1 binary format: a prefix of .bed/.bim/.fam files
  BFILE=""

  function error_exit() {
    echo "${1:-"Unknown Error"}" 1>&2
    exit 1
  }

# Input checks
  if [[ -z "${STUDY_NAME}" ]]; then
    error_exit "Please specify \$STUDY_NAME."
  fi

  if [[ -z "${ANALYST_LAST_NAME}" ]]; then
    error_exit "Please specify \$ANALYST_LAST_NAME."
  fi

  if [[ -z "${PCA_LOADINGS}" ]]; then
    error_exit "Please specify \$PCA_LOADINGS."
  fi

  if [[ -z "${PCA_AF}" ]]; then
    error_exit "Please specify \$PCA_AF."
  fi

  if [[ -n "${PFILE}" ]]; then
    input_command="--pfile ${PFILE}"
  elif [[ -n "${BFILE}" ]]; then
    input_command="--bfile ${BFILE}"
  else
    error_exit "Either \$PFILE or \$BFILE should be specified"
  fi

# Run plink2 --score
  plink2 \
    ${input_command} \
    --score ${PCA_LOADINGS} \
            variance-standardize \
            cols=-scoreavgs,+scoresums \
            list-variants \
            header-read \
    --score-col-nums 3-12 \
    --read-freq ${PCA_AF} \
    --out ${OUTNAME}

# The score file does not have FID (=0)
  awk '
  {
    if (NR==1) $1="#FID IID"
    else $1=0" "$1
    print
  }' ${OUTNAME}.sscore > ${dir}/snpid.sscore
  ln -sf ${OUTNAME}.sscore.vars ${dir}/snpid.sscore.vars

# Our PCs are PC_# rather than PC#
# The phenotype and covariate file missed FID (=0) column
  awk '
  {
    if (NR==1)
    {
      $1="FID I"$1
      gsub(/PC_/,"PC",$0)
    }
    else $1=0" " $1
  };1' ${dir}/work/INTERVAL-covid.txt | \
  tr ' ' '\t' > ${dir}/snpid.txt
  cut -f1-3 ${dir}/snpid.txt > ${dir}/snpid.pheno
  cut -f3 --complement ${dir}/snpid.txt > ${dir}/snpid.covars

  stata -b do ethnic.do

  Rscript ${PCA_projection}/plot_projected_pc.R \
        --sscore ${prefix}/work/snpid.sscore \
        --phenotype-file ${prefix}/work/snpid.pheno \
        --phenotype-col SARS_CoV \
        --covariate-file ${prefix}/work/snpid.covars \
        --pc-prefix PC \
        --pc-num 20 \
        --ancestry-file ethnic.txt \
        --ancestry-col ethnic \
        --study ${STUDY_NAME} \
        --out ${OUTNAME}
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
