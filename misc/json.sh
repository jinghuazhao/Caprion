#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export analysis=~/Caprion/analysis
export pre_qc_data=/rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics
export suffix=_dr

function tojson()
{
  export flanking=250000
  export r=$(echo ${chr}:$(expr ${pos} - ${flanking})-$(expr ${pos} + ${flanking}))
  (
    echo chromosome position variant rsid ref_allele alt_allele alt_allele_freq log_pvalue beta se
    tabix ~/Caprion/analysis/METAL_dr/${prot}_dr-1.tbl.gz ${r} | \
    awk '{
          $5=toupper($5);$4=toupper($4)
          print $1,$2,$1":"$2"_"$5"/"$4,$3,$5,$4,$6,-$12,$10,$11
    }' | \
    sort -k1,1n -k2,2n
  ) | \
  tr ' ' '\t' | \
  bgzip -f > ${prot}-${rsid}.gz
  tabix -S1 -s1 -b2 -e2 ${prot}-${rsid}.gz
  Rscript -e '
    library(dplyr)
    library(jsonlite)
    prot <- Sys.getenv("prot")
    pqtl <- Sys.getenv("rsid")
    d <- read.delim(paste0(prot,"-",pqtl,".gz"))
    sink(paste0(prot,"-",pqtl,".json"))
    toJSON(list(data=d,analysis=paste0(prot,"-",pqtl)),auto_unbox=TRUE,pretty=FALSE)
    sink()
  '
}
export -f tojson

function lz_json()
{
  awk 'NR>1{print $1,$7}' ${analysis}/work/caprion${suffix}.signals | \
  parallel -C' ' '
    export prot={1}
    export rsid={2}
    export chrpos=$(grep -w ${rsid} ~/INF/work/INTERVAL.rsid | awk "{gsub(/chr|_[A-Z]*/,\"\",\$1);print \$1}")
    export chr=$(echo ${chrpos} | cut -d ":" -f1)
    export pos=$(echo ${chrpos} | cut -d ":" -f2)
    echo ${prot} ${rsid} ${chr} ${pos}
    tojson
  '
}

function json()
{
  Rscript -e '
    library(dplyr)
    library(jsonlite)
    suffix <- Sys.getenv("suffix")
    merged_data <- list()
    flist <- dir(pattern="json")
    for (i in 1:length(flist))
    {
      s <- unlist(strsplit(flist[i],"-|[.]"))
      f <- list(ppid=paste0(s[1],"-",s[2]),data=fromJSON(flist[i])$data[c("variant","position","ref_allele","alt_allele_freq","beta","log_pvalue")])
      merged_data <- c(merged_data,list(f))
    }
    merged_json <- toJSON(merged_data)
    sink(paste0("caprion",suffix,".js"))
    cat("input=")
    writeLines(merged_json)
    sink()
  '
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


