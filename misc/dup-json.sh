#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export analysis=~/Caprion/analysis
export root=${analysis}/dup

module load ceuadmin/htslib

function gz()
{
  export flanking=250000
  export r=$(echo ${chr}:$(expr ${pos} - ${flanking})-$(expr ${pos} + ${flanking}))
  echo ${r}
  (
    echo chromosome position variant rsid ref_allele alt_allele alt_allele_freq log_pvalue beta se
    cat <(tabix ${root}/ZWK-1-${isotope}.fastGWA.gz ${r}) \
        <(tabix ${root}/ZWK-1-${isotope}-chrX.fastGWA.gz ${r}) | \
    awk '{$2="chr"$2};1' | \
    sort -k2,2 | \
    join -22 <(awk -vchr=${chr} '$1 ~ "^chr"chr":"' ${INF}/work/INTERVAL.rsid) - | \
    awk '{
          $6=toupper($6);$5=toupper($5)
          print $3,$4,$3":"$4"_"$6"/"$5,$2,$6,$5,$8,-log($11)/log(10),$9,$10
    }' | \
    sort -k1,1n -k2,2n
  ) | \
  tr ' ' '\t' | \
  bgzip -f > ${root}/json/${isotope}-${snpid}.gz
  tabix -S1 -s1 -b2 -e2 ${root}/json/${isotope}-${snpid}.gz
}

function json()
{
  Rscript -e '
    suppressMessages(library(dplyr))
    suppressMessages(library(jsonlite))
    root <- Sys.getenv("root")
    isotope <- Sys.getenv("isotope")
    pqtl <- Sys.getenv("snpid")
    d <- read.delim(file.path(root,"json",paste0(isotope,"-",pqtl,".gz")))
    f <- gzfile(file.path(root,"json","gz",paste0(isotope,"-",pqtl,".json.gz")))
    sink(f)
    toJSON(list(ppid=paste0(isotope,"-",pqtl),data=d),auto_unbox=TRUE,pretty=FALSE)
    sink()
  '
}
export -f gz
export -f json


function gz_json()
{
  awk 'NR>1{print $5,$6}' ${root}/dup.signals | \
  parallel -C' ' '
    export isotope={1}
    export snpid={2}
    IFS=\_ read chrpos a1 a2 <<<${snpid}
    IFS=\: read chr pos <<<${chrpos}
    gz
    json
    echo ${isotope} ${snpid} ${chrpos}
  '
}

function lz_json()
{
  awk 'NR>1{print $5,$6}' ${root}/dup.signals | \
  sort -k1,1n | \
  parallel -C' ' '
    export isotope={1}
    export snpid={2}
    IFS=\_ read chrpos a1 a2 <<<${snpid}
    IFS=\: read chr pos <<<${chrpos}
    echo ${isotope} ${snpid} ${chrpos}
  ' | \
  awk '!/X/' | \
  Rscript -e '
    suppressMessages(library(jsonlite))
    root <- Sys.getenv("root")
    d <- read.table("stdin",col.names=c("isotope","snpid","chrpos"))
    writeLines(toJSON(d,dataframe="values",pretty=TRUE))
    write.table(d[,-4],file=file.path(root,"json","top_hits.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  ' | gzip -f > ${root}/json/gz/top_hits.json.gz
}

gz_json
lz_json
