#!/usr/bin/bash

export caprion=~/Caprion/analysis
export work=~/Caprion/analysis/work
export TMPDIR=${HPC_WORK}/work

function setup()
{
  if [ ! -d ${caprion}/METAL/sentinels ]; then mkdir -p ${caprion}/METAL/sentinels; fi
}

function signals()
(
  cat ${caprion}/METAL/sentinels/*signals | \
  head -1 | \
  awk -v FS="\t" '{print "prot",$0}'
  cat ~/Caprion/pilot/work/caprion.varlist | \
  parallel -C' ' '
    if [ -f ${caprion}/METAL/sentinels/{}.signals ]; then
       awk -v FS="\t" -v prot={} "NR>1 {print prot,\$0}" ${caprion}/METAL/sentinels/{}.signals
    fi
    if [ -f ${caprion}/METAL/sentinels/{}-chrX.signals ]; then
       awk -v FS="\t" -v prot={} "NR>1 {print prot,\$0}" ${caprion}/METAL/sentinels/{}-chrX.signals
    fi
  '
) > ${caprion}/work/caprion.signals


function cistrans()
{
  cat <(gunzip -c ${caprion}/METAL/*-1.tbl.gz | head -1 | paste <(echo prot) -) \
      <(sed '1d;s/\t/ /g' ${work}/caprion.signals | grep -v 'X:' | \
        parallel -C' ' -j20 'zgrep -w {7} ${analysis}/METAL/{1}-1.tbl.gz | paste <(echo {1}) -') \
      <(sed '1d;s/\t/ /g' ${work}/caprion.signals | grep 'X:' | \
        parallel -C' ' -j20 'zgrep -w {7} ${analysis}/METAL/{1}-chrX-1.tbl.gz | paste <(echo {1}) -') \
      > ${work}/caprion.merge
  cut -f11,14 work/caprion.merge | sed '1d' | awk -vOFS="\t" '{printf $2" "; if($1<0) print "-"; else print "+"}' | \
  sort -k1,1 -k2,2 | uniq -c | awk -vOFS="\t" '{print $1,$2,$3}'> ${work}/caprion.dir
  Rscript -e '
  # Directions
    work <- Sys.getenv("work")
    caprion.dir <- within(read.table("work/caprion.dir",col.names=c("Count","Direction","Final")),{Direction=gsub(""," ",Direction)})
    knitr::kable(caprion.dir)
  # cis/trans classification
    signals <- read.table(file.path("work","caprion.signals"),header=TRUE)
    merged <- read.delim(file.path("work","caprion.merge"))
    names(merged)[1:4] <- c("prot","Chr","bp","SNP")
  # glist-hg19
    INF <- Sys.getenv("INF")
    glist_hg19 <- read.table(file.path(INF,"csd3","glist-hg19"),col.names=c("chr","start","end","gene")) %>%
                  filter(! chr %in% c("XY","Y"))
    X <- with(glist_hg19,chr=="X")
    glist_hg19[X,"chr"] <- "23"
    suppressMessages(library(dplyr))
    ucsc <- transmute(pQTLtools::hg19Tables,chr=gsub("chr","",X.chrom),start=chromStart,end=chromEnd,Gene=hgncSym) %>%
            select(Gene,chr,start,end)
    X <- with(ucsc,chr=="X")
    ucsc[X,"chr"] <- "23"
    added_gene <- subset(glist_hg19,gene %in%c("APOC2","APOC4","CZIB","FYB1","SEPTIN2")) %>%
                  rename(Gene=gene)
    ucsc <- bind_rows(ucsc,added_gene)
    caprion <- select(pQTLtools::caprion,Protein,Accession,Gene) %>%
               mutate(Protein=gsub("_HUMAN","",Protein)) %>%
               rename(prot=Protein)
    caprion[c(55,237,433,435),]
    subset(glist_hg19,grepl("^AMY1",gene))
    subset(glist_hg19,grepl("^C4B",gene)&chr=="6")

    pqtls <- select(merged,prot,SNP,log.P.) %>%
             rename(log10P=log.P.) %>%
             left_join(caprion) %>%
             select(Gene,SNP,log10P)
    posSNP <- select(merged,SNP,Chr,bp)
    suppressMessages(library(iBMQ))
    cis.vs.trans <- eqtlClassifier(pqtls,posSNP,ucsc,1e6)
    table(cis.vs.trans$Type)
  ##eqtlClassifier
    posGene <- select(glist_hg19,gene,chr,start,end)
    cvt <- eqtlClassifier(pqtls,posSNP,posGene,1e6)
    table(cvt$Type)
  ##cis.vs.trans.classification
    panel <- left_join(caprion,glist_hg19,by=c('Gene'='gene')) %>%
             filter(prot %in% unique(merged$prot))
    cistrans <- gap::cis.vs.trans.classification(merged,panel,"prot")
    cvt <- with(cistrans,data)
  '
}

# Number of proteins which contribute pQTLs
sed '1d' work/caprion.signals |cut -d' ' -f1 | sort -k1,1 | uniq | wc -l
# Number of unique variants
awk 'NR>1{print $7}' work/caprion.signals | sort | uniq | wc -l
# ordering of pQTLs
sed '1d;s/\t/ /g' work/caprion.signals | sort -k5,5gr
