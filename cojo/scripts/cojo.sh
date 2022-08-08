#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export bfile=${1}
export p=${2}
export r=${3}
export pr=${p}-${r}
export chr=${4}
export pos=${5}
export flanking=1e6
export start=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{start=pos-flanking;if(start<0) start=0;print start}')
export end=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{print pos+flanking}')

plink2 --bfile ${bfile} \
       --chr ${chr} --from-bp ${start} --to-bp ${end} \
       --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.8 --out results/${pr}
if [ $(grep -w ${r} results/${pr}.prune.in | wc -l) -eq 0 ]; then
   export i=$(grep -w -f results/${pr}.prune.in {bfile}.bim | \
              awk -vpos=${pos} 'function abs(x) {if (x<0) return -x; else return x;} {d=abs($4-pos);print $1, $2, $4, d}' | \
              sort -r -k4,4n | \
              awk 'NR==1 {print $2}' \
            )
   sed -i 's/'"$i"'/'"$r"'/g' results/${pr}.prune.in
fi
if [ ${chr} -eq 19 ]; then
   sort results/${pr}.prune.in | \
   join -v1 - ${INF}/work/NLRP2 > results/${pr}.prune
else
   sort results/${pr}.prune.in > results/${pr}.prune
fi
rm results/${pr}.prune.in results/${pr}.prune.out
