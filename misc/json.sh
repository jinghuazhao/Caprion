#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export analysis=~/Caprion/analysis
export pre_qc_data=/rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics
export suffix=_dr

module load ceuadmin/htslib

function gz()
{
  export flanking=250000
  export r=$(echo ${chr}:$(expr ${pos} - ${flanking})-$(expr ${pos} + ${flanking}))
  echo ${r}
  (
    echo chromosome position variant rsid ref_allele alt_allele alt_allele_freq log_pvalue beta se
    tabix ~/Caprion/analysis/METAL_dr/${prot}_dr-1.tbl.gz ${r} | \
    awk '{$3="chr"$3};1' | \
    sort -k3,3 | \
    join -23 <(awk -vchr=${chr} '$1 ~ "^chr"chr":"' ${INF}/work/INTERVAL.rsid) - | \
    awk '{
          $6=toupper($6);$5=toupper($5)
          print $3,$4,$3":"$4"_"$6"/"$5,$2,$6,$5,$7,-$13,$11,$12
    }' | \
    sort -k1,1n -k2,2n
  ) | \
  tr ' ' '\t' | \
  bgzip -f > ${analysis}/json/${prot}-${snpid}.gz
  tabix -S1 -s1 -b2 -e2 ${analysis}/json/${prot}-${snpid}.gz
}

function json()
{
  Rscript -e '
    suppressMessages(library(dplyr))
    suppressMessages(library(jsonlite))
    analysis <- Sys.getenv("analysis")
    prot <- Sys.getenv("prot")
    pqtl <- Sys.getenv("snpid")
    d <- read.delim(file.path(analysis,"json",paste0(prot,"-",pqtl,".gz")))
    sink(file.path(analysis,"json",paste0(prot,"-",pqtl,".json")))
    toJSON(list(ppid=paste0(prot,"-",pqtl),data=d),auto_unbox=TRUE,pretty=FALSE)
    sink()
  '
}
export -f gz
export -f json

function lz_json()
{
  awk 'NR>1{print $1,$7}' ${analysis}/work/caprion${suffix}.signals | \
  parallel -C' ' '
    export prot={1}
    export snpid={2}
    IFS=\_ read chrpos a1 a2 <<<${snpid}
    IFS=\: read chr pos <<<${chrpos}
    echo ${prot} ${snpid} ${chrpos}
  # gz
  # json
  ' | \
  awk '!/X/' | \
  Rscript -e '
    suppressMessages(library(jsonlite))
    analysis <- Sys.getenv("analysis")
    d <- read.table("stdin",col.names=c("prot","snpid","chrpos"))
    writeLines(toJSON(d,dataframe="values",pretty=TRUE))
    write.table(d[,-4],file=file.path(analysis,"json","top_hits.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  ' | gzip -f > ${analysis}/json/gz/top_hits.json.gz
}

function umich()
{
  curl https://portaldev.sph.umich.edu/api/v1/statistic/single/ > single.json
  cat single.json | jq
  curl -G "https://portaldev.sph.umich.edu/api/v1/statistic/phewas/?build=GRCh37&format=objects" \
       --data-urlencode "filter=variant eq '10:114758349_C/T'" > t2d.json
  curl -G "https://portaldev.sph.umich.edu/api/v1/statistic/single/results/" \
       --data-urlencode "filter=analysis in 45 and chromosome in '10' and position ge 114258349 and position le 115258349" \
       --data-urlencode "fields=variant, position, log_pvalue"  --data-urlencode "sort=log_pvalue" > tcf7l2.json
  Rsript -e '
   options(width=200)
   library(jsonlite)
   single <- fromJSON("single.json")
   names(single)
   data.frame(single$data)[45,]
  '
}

lz_json

function legacy()
# code for caprion_dr.js
# toJSON(list(ppid=paste0(prot,"-",pqtl),data=d,analysis=paste0(prot,"-",pqtl)),auto_unbox=TRUE,pretty=FALSE)
{
  Rscrpit -e '
    library(dplyr)
    library(jsonlite)
    analysis <- Sys.getenv("analysis")
    suffix <- Sys.getenv("suffix")
    vars <- c("variant","position","ref_allele","alt_allele_freq","beta","log_pvalue")
    merged_data <- list()
    jsonlist <- dir(file.path(analysis,"json"),pattern="json")
    flist <- grep("A1BG|ACE",jsonlist,value=TRUE)
    for (i in 1:length(flist))
    {
      s <- unlist(strsplit(flist[i],"-|[.]"))
      f <- list(ppid=paste0(s[1],"-",s[2]),data=fromJSON(file.path(analysis,"json",flist[i]))$data[vars])
      merged_data <- c(merged_data,list(f))
    }
    merged_json <- toJSON(merged_data)
    sink(file.path(analysis,"json",paste0("caprion",suffix,".js")))
    cat("input=")
    writeLines(merged_json)
    sink()
    decouple <- function(f)
    {
      d <- fromJSON(f)
      for (i in 1:nrow(d))
      {
        di <- with(d,list(ppid=ppid[[i]],data=data[[i]]))
        sink(paste0(i-1,".json"))
        writeLines(toJSON(di))
        sink()
      }
    }
  '
}
