#!/usr/bin/bash

export caprion=~/Caprion
export TMPDIR=${HPC_WORK}/work

if [ ! -d ${caprion}/analysis/METAL ]; then mkdir ${caprion}/analysis/METAL; fi

function METAL_list()
# build the complete list of files
{
  ls ${caprion}/analysis/work/pgwas/caprion-?-*.fastGWA.gz | \
  xargs -l basename -s .fastGWA.gz | \
  awk -vdir=${caprion}/analysis/work/pgwas -vOFS="\t" '{s=substr($1,9,1);p=$1;gsub("caprion-[0-9]-","",p);print s, p, dir"/"$1".fastGWA.gz"}' | \
  sort -k2,2 -k1,1n > ${caprion}/analysis/METAL/METAL.list
}

function METAL_list_pilot()
# build the complete list of files
{
  ls ${caprion}/pilot/work/caprion-?-*.fastGWA.gz | \
  xargs -l basename -s .fastGWA.gz | \
  awk -vdir=${caprion}/pilot/work -vOFS="\t" '{s=substr($1,9,1);p=$1;gsub("caprion-[0-9]-","",p);print s, p, dir"/"$1".fastGWA.gz"}' | \
  sort -k2,2 -k1,1n > ${caprion}/analysis/METAL/METAL.list
}

function METAL_files_suffix()
# generate individual METAL command files
{
  export suffix=${1}
  for p in $(cat ${caprion}/pilot/work/caprion.varlist)
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
     echo ADDFILTER AF1 ">=" 0.0001
     echo ADDFILTER AF1 "<=" 0.9999
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
     echo OUTFILE ${caprion}/analysis/METAL/$p${suffix}- .tbl
     awk '$2==token {print "PROCESS", $3}' token=${p}${suffix} ${caprion}/analysis/METAL/METAL.list
     echo ANALYZE HETEROGENEITY
     echo CLEAR
  ) > ${caprion}/analysis/METAL/${p}${suffix}.metal
  done
}

function METAL_files()
{
  METAL_files_suffix
  METAL_files_suffix -chrX
}

function METAL_analysis()
# conduct the analysis
# module load metal/2011-03-25
{
  export rt=${caprion}/analysis/METAL
  ls $rt/*.metal | \
  xargs -l basename -s .metal | \
  parallel -j1 --env rt -C' ' '
    metal ${rt}/{}.metal 2>&1 | \
    tee ${rt}/{}-1.tbl.log;
    gzip -f ${rt}/{}-1.tbl
  '
}

$1
