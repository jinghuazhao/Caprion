# 25-1-2021 JHZ

export caprion=${HOME}/Caprion
qctool -g ${caprion}/data/caprion-#.bgen -s SomaLogic.sample -ofiletype binary_ped -og ${caprion}/data/caprion.bgen
plink --bfile ${caprion}/data/caprion.bgen --maf 0.01 --make-bed --out ${caprion}/data/caprion.01
awk '{$1=$2};1' ${caprion}/data/caprion.01.fam > ${caprion}/data/caprion.fam
cut -f2 ${caprion}/data/caprion.01.bim > ${caprion}/data/caprion.01.snpids
qctool -g ${caprion}/data/caprion-#.bgen -s SomaLogic.sample -og ${caprion}/data/caprion.01.bgen -bgen-bits 8 -incl-snpids ${caprion}/data/caprion.01.snpids
bgenix -g ${caprion}/data/caprion.01.bgen -index -clobber
