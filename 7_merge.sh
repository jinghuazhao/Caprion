#!/usr/bin/bash

export caprion=~/Caprion/analysis
export TMPDIR=${HPC_WORK}/work

function setup()
{
  if [ ! -d ${caprion}/METAL/sentinels ]; then mkdir -p ${caprion}/METAL/sentinels; fi
}

function signals()
(
  cat ${caprion}/METAL/sentinels/*signals | \
  head -1 | \
  awk -v FS="\t" '{print "prot",$0}'
  cat ~/Caprion/pilot/work/caprion.varlist | \
  parallel -C' ' '
    if [ -f ${caprion}/METAL/sentinels/{}.signals ]; then
       awk -v FS="\t" -v prot={} "NR>1 {print prot,\$0}" ${caprion}/METAL/sentinels/{}.signals
    fi
    if [ -f ${caprion}/METAL/sentinels/{}-chrX.signals ]; then
       awk -v FS="\t" -v prot={} "NR>1 {print prot,\$0}" ${caprion}/METAL/sentinels/{}-chrX.signals
    fi
  '
) > ${caprion}/work/caprion.signals
