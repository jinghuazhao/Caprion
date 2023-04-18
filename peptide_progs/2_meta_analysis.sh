#!/usr/bin/bash

export root=~/Caprion/analysis/peptide/${protein}
export TMPDIR=${HPC_WORK}/work

if [ ! -d ${root}/METAL ]; then mkdir ${root}/METAL; fi

function METAL_list()
# build the complete list of files
{
  ls ${root}/${protein}-?-*.fastGWA.gz | \
  xargs -l basename -s .fastGWA.gz | \
  tr '-' '\t' | \
  awk -vdir=${root} -vOFS="\t" '
  {
     protein=$1;
     batch=$2;
     isotope=$3
     chrx=$4
     print protein, isotope, batch, chrx, dir"/"protein"-"batch"-"isotope".fastGWA.gz"
  }' | \
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
     echo OUTFILE ${root}/METAL/${isotope}-${suffix}- .tbl
     awk '$2==isotope && $4==suffix {print "PROCESS", $5}' FS="\t" isotope=${isotope} suffix=${suffix} ${root}/METAL/METAL.list
     echo ANALYZE HETEROGENEITY
     echo CLEAR
  ) > ${root}/METAL/${isotope}-${suffix}.metal
  done
}

function METAL_files()
{
  METAL_files_suffix
  METAL_files_suffix chrX
}

function METAL_analysis()
{
  export rt=${root}/METAL
  ls $rt/*.metal | \
  xargs -l basename -s .metal | \
  parallel -j1 --env rt -C' ' '
    metal ${rt}/{}.metal 2>&1 | \
    tee ${rt}/{}-1.tbl.log;
    gzip -f ${rt}/{}-1.tbl
  '
}

# A specific function name
$1
