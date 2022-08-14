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

function gcta()
{
  cat <(echo SNP A1 A2 freq b se p N) \
      <(gunzip -c METAL/sentinels/${p}.p.gz | awk '{print $3,toupper($4),toupper($5),$6,$10,$11,10^$12,$18}') > work/${p}.ma
  cut -d' ' -f1 work/${p}.ma | sed '1d' > work/${p}.rsid
  plink2 --bgen data/chr${chr}.bgen ref-last --sample data/caprion.sample \
         --extract work/${p}.rsid --export ind-major-bed --maf 0.01 --out work/chr${chr}
  gcta-1.9 --bfile work/chr${chr} \
           --cojo-file work/${p}.ma \
           --cojo-slct \
           --cojo-p 5e-8 \
           --cojo-collinear 0.9 \
           --out results/${pr}
}

gcta

function dosage()
{
  plink2 --bfile ${bfile} \
         --chr ${chr} --from-bp ${start} --to-bp ${end} \
         --geno 0.1 --mind 0.1 --maf 0.01 --indep-pairwise 1000kb 1 0.1 --out results/${pr}
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
  plink-2 --bgen data/chr${chr}.bgen ref-last \
          --sample data/caprion.sample \
          --extract prune/${pr}.prune \
          --recode A include-alt \
          --out results/${pr}.dosage
}

# mkdir data
# cd data
# cat ~/Caprion/pilot/work/caprion.bgenlist | xargs -l -I {} ln -s {}
# ln -s ~/Caprion/pilot/work/caprion.sample
# cd -

# legacy

function regress()
{
  Rscript -e '
    p <- Sys.getenv("p")
    pr <- Sys.getenv("pr")
    rsids <- scan(pr,"")
    chunks <- split(rsids, ceiling(seq_along(rsids)/snakemake@params[["chunksize"]]))
    load("reports/edata_batch.rda")
    # https://www.statology.org/stepwise-regression-r/
    lapply(chunks, {
      intercept_only <- lm(p ~ 1, data=edata)
      all <- lm(p ~ ., data=edata)
      both <- step(intercept_only, direction='both', scope=formula(all), trace=0)
      both$anova
      summary(both)
    })
  '
}
