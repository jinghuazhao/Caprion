#!/usr/bin/bash

function METAL_list()
# build the complete list of files
{
  ls ${root}/${protein}-?-*.fastGWA.gz | \
  xargs -l basename -s .fastGWA.gz | \
  tr '-' '\t' | \
  awk -vdir=${root} '
  {
     protein=$1;
     batch=$2;
     isotope=$3
     chrX=$4
     printf "%s,%d,%d,%s,%s", protein, isotope, batch, chrX, dir"/"protein"-"batch"-"isotope;
     if (chrX=="") print ".fastGWA.gz"; else print "-chrX.fastGWA.gz";
  }' | \
  tr ',' '\t' | \
  sort -t$'\t' -k2,2n -k4,4 -k3,3n > ${root}/METAL/METAL.list
}

function METAL_files()
# generate individual METAL command files
{
  for isotope in $(head -1 ${root}/${protein}.pheno | cut -d' ' -f1-2 --complement)
  do
  (
     echo SEPARATOR TAB
     echo COLUMNCOUNTING STRICT
     echo CHROMOSOMELABEL CHR
     echo POSITIONLABEL POS
     echo CUSTOMVARIABLE N
     echo LABEL N as N
     echo TRACKPOSITIONS ON
     echo AVERAGEFREQ ON
     echo MINMAXFREQ ON
     echo ADDFILTER AF1 ">=" 0.01
     echo ADDFILTER AF1 "<=" 0.99
     echo MARKERLABEL SNP
     echo ALLELELABELS A1 A2
     echo EFFECTLABEL BETA
     echo PVALUELABEL P
     echo WEIGHTLABEL N
     echo FREQLABEL AF1
     echo STDERRLABEL SE
     echo SCHEME STDERR
     echo EFFECT_PRINT_PRECISION 8
     echo STDERR_PRINT_PRECISION 8
     echo GENOMICCONTROL OFF
     echo LOGPVALUE ON
     echo OUTFILE ${root}/METAL/${isotope}- .tbl
     awk '$2==isotope && ($4=="") {print "PROCESS", $5}' FS="\t" isotope=${isotope} ${root}/METAL/METAL.list
     awk '$2==isotope && ($4=="chrX") {print "PROCESS", $5}' FS="\t" isotope=${isotope} ${root}/METAL/METAL.list
     echo ANALYZE HETEROGENEITY
     echo CLEAR
  ) > ${root}/METAL/${isotope}.metal
  done
}

function METAL_analysis_parallel()
{
  export rt=${root}/METAL
  ls $rt/*.metal | \
  xargs -l basename -s .metal | \
  parallel -j3 --env rt -C' ' '
    metal ${rt}/{}.metal 2>&1 | \
    tee ${rt}/{}-1.tbl.log;
    gzip -f ${rt}/{}-1.tbl
  '
}

function METAL_analysis_sbatch()
{
cat << 'EOL' > ${root}/${protein}-METAL.sb
#!/usr/bin/bash

#SBATCH --job-name=_MA-LABEL
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition cclake-himem
#SBATCH --array=1-RUNS
#SBATCH --output=PROTEIN-METAL_%A_%a.o
#SBATCH --error=PROTEIN-METAL_%A_%a.e
#SBATCH --export ALL

export TMPDIR=${HPC_WORK}/work
export rt=ROOT/METAL

function metal_parallel()
{
  ls $rt/*.metal | \
  xargs -l basename -s .metal | \
  parallel -j3 --env rt -C' ' '
      metal ${rt}/{}.metal 2>&1 | \
      tee ${rt}/{}-1.tbl.log;
      bgzip -f ${rt}/{}-1.tbl
  '
}

function metal_array()
{
  ls $rt/*.metal | \
  xargs -l basename -s .metal | \
  awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]' | \
  parallel -j3 --env rt -C' ' '
      metal ${rt}/{}.metal 2>&1 | \
      tee ${rt}/{}-1.tbl.log;
      bgzip -f ${rt}/{}-1.tbl
  '
}

metal_array

EOL
sed -i "s|PROTEIN|${root}/${protein}|;s|LABEL|${protein}|;s|ROOT|${root}|;s|RUNS|${N}|" ${root}/${protein}-METAL.sb
sbatch ${root}/${protein}-METAL.sb
#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=cclake
#metal ${rt}/${isotope}.metal 2>&1 | tee ${rt}/${isotope}-1.tbl.log;gzip -f ${rt}/${isotope}-1.tbl
#metal ${rt}/${isotope}-chrX.metal 2>&1 | tee ${rt}/${isotope}-chrX-1.tbl.log; gzip -f ${rt}/${isotope}-chrX-1.tbl
}

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export signals=${analysis}/work/caprion${suffix}.signals

# only those with pQTLs
export n_with_signals=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | wc -l)
for i in $(seq ${n_with_signals})
do
  export signal_index=${i}
  export protein=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | awk 'NR==ENVIRON["signal_index"]')
  export root=~/Caprion/analysis/peptide/${protein}
  export pheno=${analysis}/peptide/${protein}/${protein}.pheno
  export N=$(awk 'NR==1{print NF-2}' ${pheno})
  echo ${signal_index}, ${protein}
  if [ ! -d ${root}/METAL ]; then mkdir ${root}/METAL; fi
  METAL_list
  METAL_files
  METAL_analysis_sbatch
done
