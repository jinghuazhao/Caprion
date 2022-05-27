#!/usr/bin/bash

export caprion=~/Caprion/analysis
export work=~/Caprion/pilot/work
export pilot=~/Caprion/pilot/data
export hatk=~/hpc-work/HATK/

function extract_hla()
# Region as used in the SCALLOP-INF project along with a toy data
{
  plink --bfile ${pilot}/merged_imputation --chr 6 --from-bp 25392021 --to-bp 33392022 \
        --keep ${work}/caprion.id2 \
        --make-bed --out ${caprion}/work/hla
  cat <(echo FID IID sex) \
      <(awk '{print $1,$2,$5-1}' ${caprion}/work/hla.fam) > ${caprion}/work/hla.pheno
}

function hla2hped()
{
  Rscript -e '
    caprion <- Sys.getenv("caprion")
    for (gene in c("A","B","C","DPB1","DQA1","DQB1","DRB1"))
    {
        load(file.path(caprion,"HLA","HIBAG","imputed_results",paste0("hla-",gene,".rda")))
        out <- pred$value[c("sample.id","allele1","allele2","prob")]
        write.table(out,file.path(caprion,"HLA","HIBAG",paste0("hla-",gene,".txt")),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
    }
  '
  source ~/COVID-19/py37/bin/activate
  python ${hatk}/HATK.py \
         --hla2hped \
         --platform HIBAG \
         --out ${caprion}/HLA/HIBAG/hla_HIBAG \
         --rhped \
           ${caprion}/HLA/HIBAG/hla-A.txt \
           ${caprion}/HLA/HIBAG/hla-B.txt \
           ${caprion}/HLA/HIBAG/hla-C.txt \
           NA \
           ${caprion}/HLA/HIBAG/hla-DPB1.txt \
           ${caprion}/HLA/HIBAG/hla-DQA1.txt \
           ${caprion}/HLA/HIBAG/hla-DQB1.txt \
           ${caprion}/HLA/HIBAG/hla-DRB1.txt
}

# HLA imputation:
# SNP2HLA, CookHLA, HIBAG are now all part of the following SLURM script,
# sbatch 8_hla.sb
# HIBAG is specifically run through from the following R script,
# R --no-save < 8_hla.R
# CookHLA also uses SNP2HLA utility for making reference.
