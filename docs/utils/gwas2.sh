#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot

. /etc/profile.d/modules.sh 
module load ceuadmin/stata

for trait in RCN3_442625488_VADQDGDSMATR_invn RCN3_442666668_EVAKEFDQLTPEESQAR_invn RCN3_All_invn RCN3_DR_invn
do
  export trait=${trait}
  if [ ! -d ${caprion}/${trait} ]; then mkdir ${caprion}/${trait}; fi
  stata-mp -b do utils/gwas2.do
done
