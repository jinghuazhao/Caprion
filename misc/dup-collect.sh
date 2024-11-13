#!/usr/bin/bash

function step2_pqtl_collect()
{
cat <<'EOL'> ${root}/dup-step2.sb
#!/usr/bin/bash

#SBATCH --job-name=_dup2
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem

#SBATCH --export ALL
#SBATCH --output=ROOT/sentinels/slurm/_step2_%A_%a.o
#SBATCH --error=ROOT/sentinels/slurm/_step2_%A_%a.e

export TMPDIR=${HPC_WORK}/work

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/R
module load samtools/1.13/gcc/zwxn7ug3
module load perl/5.26.3_system/gcc-8.4.1-4cl2czq
module load libiconv/1.16/intel/64iicvbf
module load ceuadmin/ensembl-vep/111-icelake

function pqtls()
(
  cat ${root}/sentinels/*signals | \
  head -1 | \
  awk -v FS="\t" '{print "isotope",$0}'
  head -1 ${root}/ZWK.pheno | cut -d' ' -f1-2 --complement | tr ' ' '\n' | \
  parallel -C' ' '
    if [ -f ${root}/sentinels/{}.signals ]; then
       awk -v FS="\t" -v OFS="\t" -v isotope={} "NR>1 {print isotope,\$0}" ${root}/sentinels/{}.signals
    fi
  '
) > ${root}/dup.signals

function merge()
{
  cat <(gunzip -c ${root}/*fastGWA.gz | head -1 | paste <(echo prot) -) \
      <(sed '1d;s/\t/ /g' ${root}/${isotope}.signals | \
        parallel -C' ' -j10 'zgrep -w {7} ${root}/ZWK-1-{1}.fastGWA.gz | paste <(echo {1}) -') \
      <(sed '1d;s/\t/ /g' ${root}/${isotope}.signals | awk '$2==23' | \
        parallel -C' ' -j10 'zgrep -w {7} ${root}/ZWK-1-{1}-chrX.fastGWA.gz | paste <(echo {1}) -') \
      > ${root}/dup.merge
  cut -f11,14 ${root}/${isotope}.merge | sed '1d' | awk -vOFS="\t" '{printf $2" "; if($1<0) print "-"; else print "+"}' | \
  sort -k1,1 -k2,2 | uniq -c | awk -vOFS="\t" '{print $1,$2,$3}'> ${root}/dup.dir
}

function vep_annotate()
{
  awk -v FS="," 'NR>1{print $5}' ${root}/dup.signals | \
  sort -k1,1 | \
  uniq | \
  parallel -C' ' '
    export isotope={}
    (
      echo "##fileformat=VCFv4.0"
      echo "#CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO"
      awk -vFS="," "\$5==ENVIRON[\"isotope\"] {print \$2}" ${cvt} | \
      sort -k1,1 | \
      zgrep -f - -w ${root}/{}-1.tbl.gz | \
      cut -f1-5 | \
      awk "{gsub(/23/,\"X\",\$1);print \$1,\$2,\$3,toupper(\$4),toupper(\$5),\".\",\".\",\".\"}"
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
      echo chromosome position nearest_gene_name cistrans
      awk -vFS="," "\$5==ENVIRON[\"isotope\"] {print \$2,\$9,\$10,\$11}" ${cvt} | \
      sort -k1,1 | \
      join - <(awk "!/#/{print \$1,\$21}" ${root}/vep/{}.tab | sort -k1,1) | \
      awk "{print \$2,\$3,\$5,\$4}" | \
      sort -k1,1n -k2,2n | \
      uniq
    ) > ${root}/vep/{}.txt
  '
}

for cmd in pqtls merge vep_annotate; do $cmd; done
EOL

sed -i "s|ROOT|${root}|" ${root}/dup-step2.sb
}

source setup.sh
step2_pqtl_collect
export root=~/Caprion/analysis/dup
sbatch ${root}/dup-step2.sb
