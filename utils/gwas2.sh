#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot

. /etc/profile.d/modules.sh 
module load ceuadmin/stata

export full=(RCN3_442625488_VADQDGDSMATR RCN3_442666668_EVAKEFDQLTPEESQAR RCN3_All RCN3_DR)
export abbrev=(RCN3_44262548~n RCN3_44266666~R RCN3_All_invn RCN3_DR_invn)

sed '2,3d' ${caprion}/data2/phase2.sample > ${caprion}/data2/phase2.dat

for i in {0..3} 
do
  export y=${full[$i]}_invn
  export trait=${abbrev[$i]}
  if [ ! -d ${caprion}/${y} ]; then mkdir ${caprion}/${y}; fi
  echo ${y} -- ${trait}
  stata-mp -b do utils/gwas2.do
done
