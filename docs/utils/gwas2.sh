#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot

. /etc/profile.d/modules.sh 
module load ceuadmin/stata

function setup()
{
  sed '2,3d' ${caprion}/data2/phase2.sample > ${caprion}/data2/phase2.dat
  cd data2
  for chr in {1..22}
  do
    ln -sf caprion-${chr}.bgen chr${chr}.bgen
#   cp caprion-${chr}.bgen.bgi chr${chr}.bgen.bgi
    ln -sf ${HPC_WORK}/data/interval/chr${chr}.bgen.bgi
  done
  ln -sf ${HPC_WORK}/data/interval/SNPinfo.dta.gz
  ln -sf ${HPC_WORK}/data/interval/Chunks_15.dta
  cut -d' ' -f1-3 phase2.sample > interval.sample
  sed '1,2d' interval.sample  > sample_info.txt
  stata <<\ \ END
    insheet id idorder missing using sample_info.txt, delim(" ")
    format id %15.0g
    format idorder %15.0g
    gzsave sample_info, replace
  END
  rm sample_info.txt
  cd -
}

export full=(RCN3_442625488_VADQDGDSMATR RCN3_442666668_EVAKEFDQLTPEESQAR RCN3_All RCN3_DR)
export abbrev=(RCN3_44262548~n RCN3_44266666~R RCN3_All_invn RCN3_DR_invn)

for i in {0..3} 
do
  export y=${full[$i]}_invn
  export trait=${abbrev[$i]}
  if [ ! -d ${caprion}/${y} ]; then mkdir ${caprion}/${y}; fi
  echo ${y} -- ${trait}
  stata-mp -b do ${caprion}/utils/gwas2.do
done
