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
        parallel -C' ' -j20 'zgrep -w {7} ${caprion}/METAL/{1}-1.tbl.gz | paste <(echo {1}) -') \
      <(sed '1d;s/\t/ /g' ${work}/caprion.signals | grep 'X:' | \
        parallel -C' ' -j20 'zgrep -w {7} ${caprion}/METAL/{1}-chrX-1.tbl.gz | paste <(echo {1}) -') \
      > ${work}/caprion.merge
  cut -f11,14 ${work}/caprion.merge | sed '1d' | awk -vOFS="\t" '{printf $2" "; if($1<0) print "-"; else print "+"}' | \
  sort -k1,1 -k2,2 | uniq -c | awk -vOFS="\t" '{print $1,$2,$3}'> ${work}/caprion.dir
}

function cistrans()
{
  Rscript -e '
    options(width=120)
    suppressMessages(library(dplyr))
    suppressMessages(library(gap))
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
    caprion_modified[390,"Gene"] <- "FAM198A"
    caprion_modified[433,"Gene"] <- "HIST"
    caprion_modified[435,"Gene"] <- "HBA"
    caprion_modified[845,"Gene"] <- "SEPT2"
    APOC <- subset(glist_hg19,gene %in%c("APOC2","APOC4")) %>%
            rename(Gene=gene)
    ucsc_modified <- bind_rows(ucsc,APOC,AMY,C4B,HIST,HBA)
    pqtls <- select(merged,prot,SNP,log.P.) %>%
             mutate(log10p=-log.P.) %>%
             left_join(caprion_modified) %>%
             select(Gene,SNP,prot,log10p)
    posSNP <- select(merged,SNP,Chr,bp)
    cis.vs.trans <- qtlClassifier(pqtls,posSNP,ucsc_modified,1e6) %>%
                    mutate(geneChrom=as.integer(geneChrom),cis=if_else(Type=="cis",TRUE,FALSE))
    table(cis.vs.trans$Type)
    write.csv(cis.vs.trans,file=file.path(work,"caprion.cis.vs.trans"),row.names=FALSE,quote=FALSE)
    png(file.path(work,"caprion.pqtl2d.png"),width=12,height=10,unit="in",res=300)
    r <- qtl2dplot(cis.vs.trans,chrlen=gap::hg19,snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                   gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                   TSS=TRUE,cis="cis",plot=TRUE,cex.labels=0.6,cex.points=0.6,
                   xlab="pQTL position",ylab="Gene position")
    dev.off()
    r <- qtl2dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                     snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                     gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                     TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                     xlab="pQTL position",ylab="Gene position")
    htmlwidgets::saveWidget(r,file=file.path(work,"caprion.pqtl2dplotly.html"))
    r <- qtl3dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                     snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                     gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                     TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                     xlab="pQTL position",ylab="Gene position")
    htmlwidgets::saveWidget(r,file=file.path(work,"caprion.pqtl3dplotly.html"))
  '
}

function mean()
{
  export caprion=~/Caprion
  awk '{gsub(/NA/,"0",$NF);print}' ${caprion}/pilot/work/caprion.sample > ${caprion}/analysis/work/caprion.sample
}

function fp()
{
  export analysis=~/Caprion/analysis
  cp  ${analysis}/work/caprion.merge ${analysis}/work/tbl.tsv
  cut -f2-4 ${analysis}/work/tbl.tsv | \
  awk 'NR>1' | \
  sort -k1,2n | \
  uniq | \
  awk -vOFS="\t" '{print $1":"$2,$3}' > ${analysis}/work/rsid.tsv
  (
    gunzip -c ${analysis}/pgwas/caprion-*fastGWA.gz | head -1
    awk 'NR>1' ${analysis}/work/tbl.tsv | \
    cut -f1,4,14 --output-delimiter=' ' | \
    parallel -j10 -C' ' '
      export direction=$(zgrep -w {2} ${analysis}/METAL/{1}-1.tbl.gz | cut -f13)
      let j=1
      for i in $(grep "Input File" ${analysis}/METAL/{1}-1.tbl.info | cut -d" " -f7)
      do
         export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
         if [ "$n" != "?" ]; then zgrep -H -w {2} $i; fi
         let j=$j+1
      done
  '
  ) | \
  sed 's|/data/jinhua/INF/sumstats||g;s/.gz//g' > ${analysis}/work/all.tsv
  Rscript -e '
    require(gap)
    require(dplyr)
    tbl <- read.delim("~/Caprion/analysis/work/tbl.tsv") %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position)) %>%
           arrange(prot,SNP)
    all <- read.delim("~/Caprion/analysis/work/all.tsv") %>%
           rename(EFFECT_ALLELE=A1,REFERENCE_ALLELE=A2) %>%
           mutate(CHR=gsub("/home/jhz22/Caprion/analysis/work/pgwas/caprion-|.fastGWA","",CHR)) %>%
           mutate(batch_prot_chr=strsplit(CHR,"-|:"),
                  batch=unlist(lapply(batch_prot_chr,"[",1)),
                  prot=unlist(lapply(batch_prot_chr,"[",2)),
                  CHR=unlist(lapply(batch_prot_chr,"[",3))) %>%
           mutate(MarkerName=paste0(CHR,":",POS),
                  study=case_when(batch == batch[1] ~ "1. ZWK",
                                  batch == batch[2] ~ "2. ZYQ",
                                  batch == batch[3] ~ "3. UDP",
                                  TRUE ~ "---")) %>%
           select(-batch_prot_chr)
    rsid <- read.table("~/Caprion/analysis/work/rsid.tsv",col.names=c("MarkerName","rsid"))
    pdf("~/Caprion/analysis/work/fp.pdf",width=10,height=8)
    METAL_forestplot(tbl,all,rsid)
    dev.off()
  '
}

function HetISq()
# Code extracted from caprion.Rmd
{
  Rscript -e '
    suppressMessages(require(dplyr))
    all <- read.delim("~/Caprion/analysis/work/all.tsv") %>%
           mutate(CHR=gsub("/home/jhz22/Caprion/analysis/work/pgwas/caprion-|.fastGWA","",CHR)) %>%
           mutate(batch_prot_chr=strsplit(CHR,"-|:"),
                  batch=unlist(lapply(batch_prot_chr,"[",1)),
                  prot=unlist(lapply(batch_prot_chr,"[",2)),
                  CHR=unlist(lapply(batch_prot_chr,"[",3))) %>%
           mutate(MarkerName=paste0(CHR,":",POS),
                  Batch=case_when(batch == batch[1] ~ "1. ZWK",
                                  batch == batch[2] ~ "2. ZYQ",
                                  batch == batch[3] ~ "3. UDP",
                                  TRUE ~ "---"),
                  direction=case_when(sign(BETA) == -1 ~ "-", sign(BETA) == 1 ~ "+", sign(BETA) == 0 ~ "0", TRUE ~ "---")) %>%
           select(Batch,prot,-batch_prot_chr,MarkerName,SNP,A1,A2,N,AF1,BETA,SE,P,INFO,direction)
    b1 <- subset(all,Batch=="1. ZWK")
    names(b1) <- paste0(names(all),".ZWK")
    b1 <- rename(b1, prot=prot.ZWK, SNP=SNP.ZWK)
    b2 <- subset(all,Batch=="2. ZYQ")
    names(b2) <- paste0(names(all),".ZYQ")
    b3 <- subset(all,Batch=="3. UDP")
    names(b3) <- paste0(names(all),".UDP")
    b <- full_join(b1,b2,by=c('prot'='prot.ZYQ','SNP'='SNP.ZYQ')) %>% full_join(b3,by=c('prot'='prot.UDP','SNP'='SNP.UDP')) %>%
         mutate(directions=gsub("NA","?",paste0(direction.ZWK,direction.ZYQ,direction.UDP))) %>%
         select(-Batch.ZWK,-Batch.ZYQ,-Batch.UDP,direction.ZWK,direction.ZYQ,direction.UDP)
    tbl <- read.delim("~/Caprion/analysis/work/tbl.tsv") %>%
           arrange(prot,MarkerName) %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position), index=1:n())
    Het <- filter(tbl,HetISq>=75) %>%
           select(prot,SNP,Direction,HetISq,index) %>%
           left_join(select(b,prot,SNP,P.ZWK,P.ZYQ,P.UDP,BETA.ZWK,BETA.ZYQ,BETA.UDP))
    write.csv(Het,file="work/HetISq75.csv",row.names=FALSE,quote=FALSE)
    write(Het[['index']],file='work/HetISq75.index',sep=",",ncolumns=nrow(Het))
  '
}

function fplz()
{
  export metal=~/Caprion/analysis/METAL
# HSPB1_rs114800762 is missing as dug by the following code.
  join -a1 <(sed '1d' work/caprion.merge | awk '{print $1"_"$4}' | sort -k1,1 ) \
           <(ls METAL/qqmanhattanlz/lz/*pdf | xargs -l basename -s .pdf | awk '{print $1,NR}') | \
  awk 'NF<2' | \
  sed 's/_/ /' | \
  parallel -C' ' 'ls METAL/qqmanhattanlz/lz/{1}*pdf'
# forest/locuszoom left-right format
  ulimit -n
  ulimit -S -n 2048
  qpdf --empty --pages $(sed '1d' ~/Caprion/analysis/work/caprion.merge | sort -k1,1 -k4,4 | cut -f1,4 --output-delimiter=' ' | \
                         parallel -C' ' 'ls $(echo METAL/qqmanhattanlz/lz/{1}_{2}.pdf | sed "s/:/_/")') -- lz2.pdf
  export npages=$(qpdf -show-npages lz2.pdf)
  qpdf --pages . 1-$npages:odd -- lz2.pdf lz.pdf
# Split files, note the naming scheme
  pdfseparate lz.pdf temp-%04d-lz.pdf
  pdfseparate ${metal}/fp/fp.pdf temp-%04d-fp.pdf
# left-right with very small file size
# Combine the final pdf
  pdfjam temp-*-*.pdf --nup 2x1 --landscape --papersize '{7in,16in}' --outfile fp+lz.pdf
  rm temp*pdf
# qpdf fp+lz.pdf --pages . $(cat HetISq75.index) -- HetISq75.pdf
  qpdf fp+lz.pdf --pages . $(sed '1d' ${caprion}.merge | sort -k1,1 -k4,4 | awk '$15>=75{printf " "NR}' | sed 's/ //;s/ /,/g') -- HetISq75.pdf
}

signals
merge
cistrans
fp

function legacy()
{
# Number of proteins which contribute pQTLs
  sed '1d' work/caprion.signals |cut -d' ' -f1 | sort -k1,1 | uniq | wc -l
# Number of unique variants
  awk 'NR>1{print $7}' work/caprion.signals | sort | uniq | wc -l
# ordering of pQTLs
  sed '1d;s/\t/ /g' work/caprion.signals | sort -k5,5gr
  Rscript -e '
  # Not quite working below.
  ##qtlClassifier
    posGene <- select(glist_hg19,gene,chr,start,end)
    cvt <- gap::qtlClassifier(pqtls,posSNP,posGene,1e6)
    table(cvt$Type)
  ##cis.vs.trans.classification
    panel <- left_join(caprion_modified,ucsc) %>%
             filter(prot %in% unique(merged$prot))
    cvt <- gap::cis.vs.trans.classification(merged,panel,"prot")
  # suppressMessages(library(iBMQ))
  # cis.vs.trans <- eqtlClassifier(pqtls,posSNP,ucsc_modified,1e6)
  '
}
