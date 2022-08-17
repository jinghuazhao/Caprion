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

function cojo()
{
  cat <(echo SNP A1 A2 freq b se p N) \
      <(gunzip -c METAL/sentinels/${p}.p.gz | awk '{print $3,toupper($4),toupper($5),$6,$10,$11,10^$12,$18}') > work/${p}.ma
  cut -d' ' -f1 work/${p}.ma | sed '1d' > work/${p}.rsid
  # single significant variant
  if [ $(wc -l work/${p}.rsid | cut -d' ' -f1) -eq 1 ]; then return 0; fi
  plink2 --bgen data/chr${chr}.bgen ref-unknown --sample data/caprion.sample \
         --extract work/${p}.rsid --geno 0.1 --mind 0.1 --indep-pairwise 1000kb 1 0.1 \
         --out results/${pr}.prune
  # only one variant left
  if [ ! -f results/${pr}.prune.prune.in ]; then return 0; fi
  if [ $(wc -l results/${pr}.prune.prune.in | cut -d' ' -f1) -eq 1 ]; then return 0; fi
  if [ $(grep -w ${r} results/${pr}.prune.prune.in | wc -l | cut -d' ' -f1) -eq 0 ]; then
     export i=$(grep -w -f results/${pr}.prune.prune.in ${bfile}.bim | \
                awk -vpos=${pos} 'function abs(x) {if (x<0) return -x; else return x;} {print $1, $2, $4, abs($4-pos)}' | \
                sort -r -k4,4n | \
                awk 'NR==1 {print $2}' \
               )
     sed -i 's/'"$i"'/'"$r"'/g' results/${pr}.prune.prune.in
  fi
  (
    if [ ${chr} -eq 19 ]; then
       sort results/${pr}.prune.prune.in | join -v1 - ${INF}/work/NLRP2
    else
       sort results/${pr}.prune.prune.in
    fi
  ) > results/${pr}.prune
# rm results/${pr}.prune.prune.in results/${pr}.prune.prune.out work/${p}.rsid
  plink2 --bgen data/chr${chr}.bgen ref-unknown \
         --sample data/caprion.sample \
         --extract results/${pr}.prune \
         --recode A include-alt \
         --out results/${pr}.dosage
  plink2 --bgen data/chr${chr}.bgen ref-unknown --sample data/caprion.sample \
         --extract results/${pr}.prune --export ind-major-bed --out work/${pr}
  gcta-1.9 --bfile work/${pr} \
           --cojo-file work/${p}.ma --extract-region-snp ${r} 1000 0.1 --maf 0.01 \
           --cojo-slct \
           --cojo-p 5e-8 \
           --cojo-collinear 0.9 \
           --out results/${pr}.gcta
  Rscript -e '
    suppressMessages(library(dplyr))
    suppressMessages(library(gap))
    suppressMessages(library(lars))
    p <- Sys.getenv("p")
    pr <- Sys.getenv("pr")
    r <- scan(paste0("results/",pr,".prune"),"")
#   chunks <- split(r, ceiling(seq_along(r)/snakemake@params[["chunksize"]]))
    raw <- read.delim(paste0("results/",pr,".dosage.raw")) %>% select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)
    load("reports/edata_batch.rda")
    load("data/pheno.rda")
    edata <- with(edata_batch,edata)
    edata[,p] <- invnormal(edata[,p])
    xp <- gsub("^([0-9])","X\\1",p)
    edata <- left_join(subset(id,!is.na(IID)),data.frame(caprion_id=rownames(edata),edata)) %>%
             select(IID,xp) %>%
             left_join(raw) %>%
             select(-IID)
    names(edata) <- gsub("^X([0-9])","\\1",names(edata))
    dosage <- select(edata,-p)
    intercept_only <- lm(edata[[p]] ~ 1, data=dosage)
    all <- lm(edata[[p]] ~ ., data=dosage)
    both <- step(intercept_only, direction="both", scope=formula(all), trace=0)
    print(both$anova)
    print(summary(both))
    m <- as.matrix(subset(edata,!is.na(edata[[p]])))
    fit <- lars(m,edata[[p]][!is.na(edata[[p]])],type="lasso")
    print(fit)
  ' > results/${pr}.lm.log
}

cojo

function setup()
{
  mkdir data
  cd data
  cat ~/Caprion/pilot/work/caprion.bgenlist | xargs -l -I {} ln -s {}
  ln -s ~/Caprion/pilot/work/caprion.sample
  ln -s /home/jhz22/Caprion/pilot/work/caprion.pheno
  Rscript -e '
  pheno <- read.delim("data/caprion.pheno")
  id <- pheno[c("IID","caprion_id")]
  save(pheno,id,file="data/pheno.rda")
  '
  cd -
# https://www.statology.org/stepwise-regression-r/
}

function checkfreq()
{
  cat results/ERAP2-rs2910686.gcta.freq.badsnps | cut -f1 | sed '1d' > badsnps
  plink2 --bfile work/ERAP2-rs2910686 --freq --missing --out badsnps
  plink2 --bfile ~/Caprion/pilot/data/caprion.01 --extract badsnps --freq --missing \
         --pheno <(awk {'if(NR==1) else $1=NR;print}' ~/Caprion/pilot/data/caprion.dat) --pheno-name Q6P179_invn --out ZWK
  plink2 --bfile ~/Caprion/pilot/data2/caprion.01 --extract badsnps --freq --missing \
         --pheno <(awk {'if(NR==1) else $1=NR;print}' ~/Caprion/pilot/data/phase2.pheno) --pheno-name ERAP2_All_invn --out ZYQ
  plink2 --bfile ~/Caprion/pilot/data3/caprion.01 --extract badsnps --freq --missing \
         --pheno <(awk {'if(NR==1) else $1=NR;print}' ~/Caprion/pilot/data/protein.pheno) --pheno-name ERAP2_invn --out ZYQ

  cat ZWK.afreq ZYQ.afreq UDP.afreq <(head -1 ZWK.afreq) <(grep -w -f badsnps badsnps.afreq) | awk -vOFS='|' '{print $1,$2,$3,$4,$5,$6}'
}
