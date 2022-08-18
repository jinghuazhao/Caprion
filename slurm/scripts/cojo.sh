#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export bfile=${1}
export p=${2}
export r=${3}
export pr=${p}-${r}
export chr=${4}
export pos=${5}

function cojo()
{
  cat <(echo SNP A1 A2 freq b se p N) \
      <(gunzip -c METAL/sentinels/${p}.p.gz | awk '{print $3,toupper($4),toupper($5),$6,$10,$11,10^$12,$18}') > work/${pr}.ma
  cut -d' ' -f1 work/${pr}.ma | sed '1d' > work/${pr}.rsid
  # singleton variant
  if [ $(wc -l work/${pr}.rsid | cut -d' ' -f1) -eq 1 ]; then return 0; fi
  plink2 --bgen data/chr${chr}.bgen ref-unknown --sample data/caprion.sample \
         --extract work/${pr}.rsid --export ind-major-bed --out work/${pr}
  gcta-1.9 --bfile work/${pr} \
           --cojo-file work/${pr}.ma --maf 0.01 --diff-freq 1 \
           --cojo-slct \
           --cojo-p 5e-8 \
           --cojo-collinear 0.9 \
           --out results/${pr}.gcta
  # for flanking window adding --extract-region-snp ${r} 1000
  if [ ! -f results/${pr}.gcta.jma.cojo ]; then return 0; fi
  sed '1d' results/${pr}.gcta.jma.cojo | cut -f2 > results/${pr}.jma
  rm work/${pr}.ma work/${pr}.rsid
  plink2 --bgen data/chr${chr}.bgen ref-unknown \
         --sample data/caprion.sample \
         --extract results/${pr}.jma \
         --recode A include-alt \
         --out results/${pr}.dosage
  Rscript -e '
    suppressMessages(library(dplyr))
    suppressMessages(library(gap))
    suppressMessages(library(lars))
    p <- Sys.getenv("p")
    pr <- Sys.getenv("pr")
    r <- scan(paste0("results/",pr,".jma"),"")
    ldr <- read.delim(paste0("results/",pr,".gcta.ldr.cojo"))
    rownames(ldr) <- ldr[["SNP"]]
    ldr <- ldr[,2:(nrow(ldr)+1)]
    raw <- read.delim(paste0("results/",pr,".dosage.raw")) %>% select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)
    load("reports/edata_batch.rda")
    load("data/pheno.rda")
    edata <- with(edata_batch,edata)
    edata[,p] <- invnormal(edata[,p])
    xp <- gsub("^([0-9])","X\\1",p)
    edata <- left_join(subset(id,!is.na(IID)),data.frame(caprion_id=rownames(edata),edata)) %>%
             select(IID,all_of(xp)) %>%
             left_join(raw) %>%
             select(-IID)
    names(edata) <- gsub("^X([0-9])","\\1",names(edata))
    dosage <- select(edata,-all_of(p))
    names(dosage) <- unlist(lapply(strsplit(names(dosage),"_"),"[",1))
    r <- gsub(":",".",unlist(lapply(strsplit(r,"_"),"[",1)))
    all <- lm(edata[[p]] ~ ., data=dosage[r])
    print(anova(all))
    print(summary(all))
    m <- as.matrix(subset(edata,!is.na(edata[[p]])))
    fit <- lars(m,edata[[p]][!is.na(edata[[p]])],type="lasso")
    print(fit)
  ' > results/${pr}.lm.log
}

cojo
