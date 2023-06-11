#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis

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

function METAL_files_suffix()
# generate individual METAL command files
{
  export suffix=${1}
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
     echo OUTFILE ${root}/METAL/${isotope}${suffix}- .tbl
     awk '$2==isotope && ("-"$4==suffix||$4==suffix) {print "PROCESS", $5}' FS="\t" isotope=${isotope} suffix=${suffix} ${root}/METAL/METAL.list
     echo ANALYZE HETEROGENEITY
     echo CLEAR
  ) > ${root}/METAL/${isotope}${suffix}.metal
  done
}

function METAL_files()
{
  METAL_files_suffix
  METAL_files_suffix -chrX
}

function METAL_analysis()
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

#SBATCH --job-name=_METAL
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio

#SBATCH --output=PROTEIN-METAL.o
#SBATCH --error=PROTEIN-METAL.e
#SBATCH --export ALL

export TMPDIR=${HPC_WORK}/work
export rt=ROOT/METAL

ls $rt/*.metal | \
xargs -l basename -s .metal | \
parallel -j3 --env rt -C' ' '
    metal ${rt}/{}.metal 2>&1 | \
    tee ${rt}/{}-1.tbl.log;
    bgzip -f ${rt}/{}-1.tbl
'

EOL
sed -i "s|PROTEIN|${root}/${protein}|;s|ROOT|${root}|" ${root}/${protein}-METAL.sb
echo sbatch ${root}/${protein}-METAL.sb
#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=cclake
#metal ${rt}/${isotope}.metal 2>&1 | tee ${rt}/${isotope}-1.tbl.log;gzip -f ${rt}/${isotope}-1.tbl
#metal ${rt}/${isotope}-chrX.metal 2>&1 | tee ${rt}/${isotope}-chrX-1.tbl.log; gzip -f ${rt}/${isotope}-chrX-1.tbl
}

for i in $(seq 987)
do
  export SLURM_ARRAY_TASK_ID=${i}
  export protein=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${pilot}/work/caprion.varlist)
  export root=~/Caprion/analysis/peptide/${protein}

  if [ ! -d ${root}/METAL ]; then mkdir ${root}/METAL; fi

# A specific function name
  METAL_list
  METAL_files
  METAL_analysis_sbatch
done
