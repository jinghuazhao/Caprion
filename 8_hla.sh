#!/usr/bin/bash

export analysis=~/Caprion/analysis
export work=~/Caprion/pilot/work
export pilot=~/Caprion/pilot/data
export hatk=~/hpc-work/HATK/
export suffix=_dr

function extract_hla()
# Region as used in the SCALLOP-INF project along with a toy data
{
  plink --bfile ${pilot}/merged_imputation --chr 6 --from-bp 25392021 --to-bp 33392022 \
        --keep ${work}/caprion.id2 \
        --make-bed --out ${analysis}/work/hla
  cat <(echo FID IID sex) \
      <(awk '{print $1,$2,$5-1}' ${analysis}/work/hla.fam) > ${analysis}/work/hla.pheno
}

function hla2hped()
{
  Rscript -e '
    analysis <- Sys.getenv("analysis")
    for (gene in c("A","B","C","DPB1","DQA1","DQB1","DRB1"))
    {
        load(file.path(analysis,"HLA","HIBAG","imputed_results",paste0("hla-",gene,".rda")))
        out <- pred$value[c("sample.id","allele1","allele2","prob")]
        write.table(out,file.path(analysis,"HLA","HIBAG",paste0("hla-",gene,".txt")),
                    col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
    }
  '
  source ~/COVID-19/py37/bin/activate
  python ${hatk}/HATK.py \
         --hla2hped \
         --platform HIBAG \
         --out ${analysis}/HLA/HIBAG/hla_HIBAG \
         --rhped \
           ${analysis}/HLA/HIBAG/hla-A.txt \
           ${analysis}/HLA/HIBAG/hla-B.txt \
           ${analysis}/HLA/HIBAG/hla-C.txt \
           NA \
           ${analysis}/HLA/HIBAG/hla-DPB1.txt \
           ${analysis}/HLA/HIBAG/hla-DQA1.txt \
           ${analysis}/HLA/HIBAG/hla-DQB1.txt \
           ${analysis}/HLA/HIBAG/hla-DRB1.txt
}

# HLA imputation:
# SNP2HLA, CookHLA, HIBAG are now all part of the following SLURM script,
# sbatch 8_hla.sb
# HIBAG is specifically run through from the following R script,
# R --no-save < 8_hla.R
# CookHLA also uses SNP2HLA utility for making reference.

function hla_signals()
{
  awk '$2==6 {print $7}' ${analysis}/work/caprion${suffix}.signals | \
  grep -f - -w ${INF}/work/INTERVAL.rsid | \
  sed 's/chr6://;s/_/\t/' | \
  awk -v start=25392021 -v end=33392022 'start <= $1 && $1<= end {print $1,$3}' | \
  sort -k2,2 | \
  join -12 -27 - <(awk '$2==6' ${analysis}/work/caprion${suffix}.signals | sort -k7,7) | \
  cut -d' ' -f3 > ${analysis}/work/hla.prot
}

# prot <- scan("hla.prot",what="")
# gene <- subset(pQTLdata::caprion[1:3],Protein %in% paste0(prot,"_HUMAN"))
# write.table(gene,file="hla.gene",quote=FALSE,row.names=FALSE,sep="\t")

function hla_tapas()
{
export CookHLA=${analysis}/HLA/CookHLA
export results=${CookHLA}/results
if [ ! -d ${results} ]; then mkdir ${results}; fi
for batch in 1 2 3
  do
  for prot in P1R7 MXRA5 COCA1 HEM2 CHIT1 KIF4A CO4A TCPQ CALU NAR3 CO4B ANXA6 HEG1 CO2 LMNA CFAB TIE1 TRY3 TENX
  do
    for exon in 2.0.5 2.1.5 2.1 3.0.5 3.1.5 3.1 4.0.5 4.1.5 4.1
    do
    python -m HLAassoc LINEAR \
           --vcf ${CookHLA}/hla_CookHLA.MHC.QC.exon${exon}.raw_imputation_out.vcf \
           --out ${results}/${prot}-${exon} \
           --pheno ${analysis}/work/caprion-${batch}${suffix}.phenotype \
           --pheno-name ${prot} \
           --hped ${CookHLA}/hla_CookHLA.MHC.HLA_IMPUTATION_OUT.hped \
           --chped ${CookHLA}/interval.imgt3320.4field.chped
    done
  done
done
}
