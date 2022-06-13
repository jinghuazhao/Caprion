#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export caprion=~/Caprion
export analysis=~/Caprion/analysis
export phenoname=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${caprion}/work/caprion.varlist)
export pilot=~/Caprion/pilot
export work=${analysis}/work

cut -d',' -f2 ${pilot}/caprion_inf1.chk | grep -f - ${INF}/work/INF1.merge
