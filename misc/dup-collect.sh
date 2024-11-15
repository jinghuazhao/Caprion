#!/usr/bin/bash

function step2_pqtl_collect()
{
cat <<'EOL'> ${root}/dup-step2.sh
export TMPDIR=${HPC_WORK}/work

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/R
module load samtools/1.13/gcc/zwxn7ug3
module load perl/5.26.3_system/gcc-8.4.1-4cl2czq
module load libiconv/1.16/intel/64iicvbf
module load ceuadmin/ensembl-vep/111-icelake

function merge()
{
  (
    cat ${root}/sentinels/*signals | head -1
    head -1 ${root}/ZWK.pheno | cut -d' ' -f1-2 --complement | tr ' ' '\n' | \
    parallel -C' ' '
      if [ -f ${root}/sentinels/{}.signals ]; then
         awk "NR>1" ${root}/sentinels/{}.signals
      fi
    '
  ) > ${root}/dup.signals

  cat <(gunzip -c ${root}/*fastGWA.gz | head -1) \
      <(sed '1d;s/\t/ /g' ${root}/dup.signals | \
        parallel -C' ' -j10 'zgrep -w {6} ${root}/ZWK-1-{5}.fastGWA.gz') \
      <(sed '1d;s/\t/ /g' ${root}/dup.signals | awk '$1=="chr23"' | \
        parallel -C' ' -j10 'zgrep -w {6} ${root}/ZWK-1-{5}-chrX.fastGWA.gz') | \
  awk 'NF>1' > ${root}/dup.merge
}

function vep_annotate()
{
  awk 'NR>1{print $5}' ${root}/dup.signals | \
  sort -k1,1 | \
  uniq | \
  parallel -C' ' '
    export isotope={}
    (
      echo "##fileformat=VCFv4.0"
      echo "#CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO"
      awk "\$5==ENVIRON[\"isotope\"] {print \$NF}" ${root}/dup.signals | \
      sort -k1,1 | \
      grep -f - -w <(zcat ${root}/ZWK-1-{}.fastGWA.gz ${root}/ZWK-1-{}-chrX.fastGWA.gz) | \
      cut -f1-5 | \
      awk "{gsub(/23/,\"X\",\$1);print \$1,\$3,\$2,toupper(\$4),toupper(\$5),\".\",\".\",\".\"}"
    ) | \
    tr " " "\t" > ${root}/vep/{}.vcf
  # VEP annotation
    cd ${HPC_WORK}/loftee
    vep --input_file ${root}/vep/{}.vcf \
        --output_file ${root}/vep/{}.tab --force_overwrite \
        --cache --dir_cache /usr/local/Cluster-Apps/ceuadmin/ensembl-vep/111-icelake/.vep \
        --offline \
        --species homo_sapiens --assembly GRCh37 --pick --nearest symbol --symbol \
        --tab
    cd -
    (
      echo chromosome position nearest_gene_name
      awk "\$5==ENVIRON[\"isotope\"] {split(\$6,chrpos,\"_\");split(chrpos[1],a,\":\");print \$6,a[1],a[2]}" ${root}/dup.signals | \
      sort -k1,1 | \
      join - <(awk "!/#/{print \$1,\$21}" ${root}/vep/{}.tab | sort -k1,1) | \
      awk "{print \$2,\$3,\$5,\$4}" | \
      sort -k1,1n -k2,2n | \
      uniq
    ) > ${root}/vep/{}.txt
  '
}

vep_annotate
EOL
}

export analysis=~/Caprion/analysis
export root=~/Caprion/analysis/dup
step2_pqtl_collect
bash ${root}/dup-step2.sh
