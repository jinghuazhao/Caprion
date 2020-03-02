# 2-3-2020 JHZ

qctool -g caprion-#.bgen -s SomaLogic.sample -ofiletype binary_ped -og caprion.bgen
plink --bfile caprion.bgen --maf 0.01 --make-bed --out caprion.01
cut -f2 caprion.01.bim > caprion.01.snpids

qctool -g caprion-#.bgen -s SomaLogic.sample -og caprion.01.bgen -bgen-bits 8 -incl-snpids caprion.01.snpids
bgenix -g caprion.01.bgen -index -clobber
