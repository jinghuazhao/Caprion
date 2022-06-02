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

function merge()
{
  cat <(gunzip -c ${caprion}/METAL/*-1.tbl.gz | head -1 | paste <(echo prot) -) \
      <(sed '1d;s/\t/ /g' ${work}/caprion.signals | grep -v 'X:' | \
        parallel -C' ' -j20 'zgrep -w {7} ${analysis}/METAL/{1}-1.tbl.gz | paste <(echo {1}) -') \
      <(sed '1d;s/\t/ /g' ${work}/caprion.signals | grep 'X:' | \
        parallel -C' ' -j20 'zgrep -w {7} ${analysis}/METAL/{1}-chrX-1.tbl.gz | paste <(echo {1}) -') \
      > ${work}/caprion.merge
  cut -f11,14 ${work}/caprion.merge | sed '1d' | awk -vOFS="\t" '{printf $2" "; if($1<0) print "-"; else print "+"}' | \
  sort -k1,1 -k2,2 | uniq -c | awk -vOFS="\t" '{print $1,$2,$3}'> ${work}/caprion.dir
}

function cistrans()
{
  Rscript -e '
    suppressMessages(library(dplyr))
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
    ucsc <- transmute(pQTLtools::hg19Tables,chr=gsub("chr","",X.chrom),start=chromStart,end=chromEnd,Gene=hgncSym) %>%
            select(Gene,chr,start,end)
    X <- with(ucsc,chr=="X")
    ucsc[X,"chr"] <- "23"
    caprion <- select(pQTLtools::caprion,Protein,Accession,Gene) %>%
               mutate(Protein=gsub("_HUMAN","",Protein)) %>%
               rename(prot=Protein)
    quadruple <- function(d,label) data.frame(Gene=label,chr=min(d$chr),start=min(d$start),end=max(d$end))
    caprion[c(55,237,433,435),]
    subset(glist_hg19,grepl("^AMY1",gene))
    subset(glist_hg19,grepl("^C4B",gene)&chr=="6")
    subset(glist_hg19,grepl("^HIST1H4|^HIST2H4[A-B]|HIST4H4",gene))
    subset(glist_hg19,grepl("^HBA",gene))
    AMY <- quadruple(subset(glist_hg19,grepl("^AMY1",gene)),label="AMY")
    C4B <- quadruple(subset(glist_hg19,grepl("^C4B",gene)&chr=="6"),label="C4B")
    HIST <- quadruple(subset(glist_hg19,grepl("^HIST1H4|^HIST2H4[A-B]|HIST4H4",gene)),label="HIST")
    HBA <- quadruple(subset(glist_hg19,grepl("^HBA",gene)),label="HBA")
    caprion_modified <- caprion
    caprion_modified[55,"Gene"] <- "AMY"
    caprion_modified[237,"Gene"] <- "C4B"
    caprion_modified[278,"Gene"] <- "C1orf123"
    caprion_modified[385,"Gene"] <- "FYB"
    caprion_modified[433,"Gene"] <- "HIST"
    caprion_modified[435,"Gene"] <- "HBA"
    caprion_modified[845,"Gene"] <- "SEPT2"
    APOC <- subset(glist_hg19,gene %in%c("APOC2","APOC4")) %>%
            rename(Gene=gene)
    ucsc_modified <- bind_rows(ucsc,APOC,AMY,C4B,HIST,HBA)
    pqtls <- select(merged,prot,SNP,log.P.) %>%
             rename(log10P=log.P.) %>%
             left_join(caprion_modified) %>%
             select(Gene,SNP,log10P)
    posSNP <- select(merged,SNP,Chr,bp)
    cis.vs.trans <- gap::qtlClassifier(pqtls[,1:2],posSNP,ucsc_modified,1e6)
    table(cis.vs.trans$Type)
  ##eqtlClassifier
    posGene <- select(glist_hg19,gene,chr,start,end)
    cvt <- eqtlClassifier(pqtls,posSNP,posGene,1e6)
    table(cvt$Type)
  ##cis.vs.trans.classification
    panel <- left_join(caprion_modified,ucsc) %>%
             filter(prot %in% unique(merged$prot))
    cvt <- gap::cis.vs.trans.classification(merged,panel,"prot")
  # suppressMessages(library(iBMQ))
  # cis.vs.trans <- eqtlClassifier(pqtls,posSNP,ucsc_modified,1e6)
  '
}

cistrans

# Number of proteins which contribute pQTLs
sed '1d' work/caprion.signals |cut -d' ' -f1 | sort -k1,1 | uniq | wc -l
# Number of unique variants
awk 'NR>1{print $7}' work/caprion.signals | sort | uniq | wc -l
# ordering of pQTLs
sed '1d;s/\t/ /g' work/caprion.signals | sort -k5,5gr
