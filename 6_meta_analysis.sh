#!/usr/bin/bash

export analysis=~/Caprion/analysis
export suffix=_dr
export TMPDIR=${HPC_WORK}/work

if [ ! -d ${analysis}/METAL${suffix} ]; then mkdir ${analysis}/METAL${suffix}; fi

function METAL_list()
# build the complete list of files
{
  ls ${analysis}/pgwas${suffix}/caprion-?-*.fastGWA.gz | \
  xargs -l basename -s .fastGWA.gz | \
  awk -vdir=${analysis}/pgwas${suffix} -vOFS="\t" '{s=substr($1,9,1);p=$1;gsub("caprion-[0-9]-","",p);print s, p, dir"/"$1".fastGWA.gz"}' | \
  sort -k2,2 -k1,1n > ${analysis}/METAL${suffix}/METAL.list
}

function METAL_files()
# generate individual METAL command files
{
  for p in $(cat ${analysis}/work/caprion${suffix}.varlist)
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
     echo ADDFILTER AF1 ">=" 0.001
     echo ADDFILTER AF1 "<=" 0.999
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
     echo OUTFILE ${analysis}/METAL${suffix}/$p${suffix}- .tbl
     awk '$2 == token {print "PROCESS", $3}' token=${p}${suffix} ${analysis}/METAL${suffix}/METAL.list
     awk '$2 == token {print "PROCESS", $3}' token=${p}${suffix}-chrX ${analysis}/METAL${suffix}/METAL.list
     echo ANALYZE HETEROGENEITY
     echo CLEAR
  ) > ${analysis}/METAL${suffix}/${p}${suffix}.metal
  done
}

function METAL_analysis()
# conduct the analysis
# module load metal/2011-03-25
{
  export rt=${analysis}/METAL${suffix}
  ls $rt/*.metal | \
  xargs -l basename -s .metal | \
  parallel -j1 --env rt -C' ' '
    metal ${rt}/{}.metal 2>&1 | \
    tee ${rt}/{}-1.tbl.log;
    cat <(head -1 ${st}/{}-1.tbl) <(sed '1d' ${rt}/{}-1.tbl | sort -k1,1n -k2,2n) | \
    bgzip -f ${rt}/{}-1.tbl.gz
    tabix -S1 -s1 -b2 -e2 -f ${rt}/{}-1.tbl.gz
    rm ${rt}/{}-1.tbl
  '
}

for cmd in METAL_list METAL_files; do $cmd; done
