#!/usr/bin/bash

export caprion=${HOME}/Caprion

function phase1()
{
  qctool -g ${caprion}/data/caprion-#.bgen -s SomaLogic.sample -ofiletype binary_ped -og ${caprion}/data/caprion.bgen
  plink --bfile ${caprion}/data/caprion.bgen --maf 0.01 --make-bed --out ${caprion}/data/caprion.01
  awk '{$1=$2};1' ${caprion}/data/caprion.01.fam > ${caprion}/data/caprion.fam
  cut -f2 ${caprion}/data/caprion.01.bim > ${caprion}/data/caprion.01.snpids
  qctool -g ${caprion}/data/caprion-#.bgen -s SomaLogic.sample -og ${caprion}/data/caprion.01.bgen -bgen-bits 8 \
         -incl-snpids ${caprion}/data/caprion.01.snpids
  bgenix -g ${caprion}/data/caprion.01.bgen -index -clobber
}

function phase2()
{
  qctool -g ${caprion}/data2/caprion-#.bgen -s ${caprion}/data2/caprion-1.samples -ofiletype binary_ped -og ${caprion}/data2/caprion.bgen
  plink --bfile ${caprion}/data2/caprion.bgen --maf 0.01 --make-bed --out ${caprion}/data2/caprion.01
  awk '{$1=$2};1' ${caprion}/data2/caprion.01.fam > ${caprion}/data2/caprion.fam
  cut -f2 ${caprion}/data2/caprion.01.bim > ${caprion}/data2/caprion.01.snpids
  qctool -g ${caprion}/data2/caprion-#.bgen -s ${caprion}/data2/caprion-1.samples -og ${caprion}/data2/caprion.01.bgen -bgen-bits 8 \
         -incl-snpids ${caprion}/data2/caprion.01.snpids
  bgenix -g ${caprion}/data2/caprion.01.bgen -index -clobber
}

phase2
