#!/usr/bin/bash

export analysis=~/Caprion/analysis
export suffix=_dr
export TMPDIR=${HPC_WORK}/work
export PERL5LIB=

module load ceuadmin/R
module load perl/5.26.3_system/gcc-8.4.1-4cl2czq

function signals()
(
  cat ${analysis}/METAL${suffix}/sentinels/*signals | \
  head -1 | \
  awk -v FS="\t" '{print "prot",$0}'
  cat ${analysis}/output/caprion${suffix}.varlist | \
  parallel -C' ' -j5 '
    if [ -f ${analysis}/METAL${suffix}/sentinels/{}${suffix}.signals ]; then
       awk -v FS="\t" -v prot={} "NR>1 {print prot,\$0}" ${analysis}/METAL${suffix}/sentinels/{}${suffix}.signals
    fi
  '
) > ${analysis}/work/caprion${suffix}.signals

function merge()
{
  cat <(cat ${analysis}/work/hdr| paste <(echo prot) -) \
      <(sed '1d;s/\t/ /g' ${analysis}/work/caprion${suffix}.signals | \
        parallel -C' ' -j20 'zgrep -w {7} ${analysis}/METAL${suffix}/{1}${suffix}-1.tbl.gz | paste <(echo {1}) -') \
      > ${analysis}/work/caprion${suffix}.merge
  cut -f11,14 ${analysis}/work/caprion${suffix}.merge | sed '1d' | awk -vOFS="\t" '{printf $2" "; if($1<0) print "-"; else print "+"}' | \
  sort -k1,1 -k2,2 | uniq -c | awk -vOFS="\t" '{print $1,$2,$3}'> ${analysis}/work/caprion${suffix}.dir
}

function cistrans()
{
  Rscript -e '
    options(width=120)
    suppressMessages(library(dplyr))
    suppressMessages(library(gap))
  # Directions
    analysis <- Sys.getenv("analysis")
    suffix=Sys.getenv("suffix")
    caprion.dir <- within(read.table(paste0(analysis,"/work/caprion",suffix,".dir"),
                                     col.names=c("Count","Direction","Final")),{Direction=gsub(""," ",Direction)})
    knitr::kable(caprion.dir)
  # cis/trans classification
    signals <- read.table(file.path(analysis,"work",paste0("caprion",suffix,".signals")),header=TRUE)
    merged <- read.delim(file.path(analysis,"work",paste0("caprion",suffix,".merge")))
    names(merged)[1:4] <- c("prot","Chr","bp","SNP")
  # glist-hg19
    INF <- Sys.getenv("INF")
    glist_hg19 <- read.table(file.path(INF,"csd3","glist-hg19"),col.names=c("chr","start","end","gene")) %>%
                  filter(! chr %in% c("XY","Y"))
    X <- with(glist_hg19,chr=="X")
    glist_hg19[X,"chr"] <- "23"
    sept7 <- filter(pQTLdata::hg19, SYMBOL == "SEPTIN7") %>%
             group_by(SYMBOL, chr) %>%
             summarize(start = min(start), end = max(end), .groups = 'drop') %>%
             transmute(Gene = as.character(SYMBOL), chr = gsub("chr", "", chr), start, end) %>%
             data.frame()
    ucsc <- transmute(pQTLdata::hg19Tables,chr=gsub("chr","",X.chrom),start=chromStart,end=chromEnd,Gene=hgncSym) %>%
            select(Gene,chr,start,end) %>%
            bind_rows(sept7)
    X <- with(ucsc,chr=="X")
    ucsc[X,"chr"] <- "23"
    caprion <- select(pQTLdata::caprion,Protein,Accession,Gene) %>%
               mutate(Protein=gsub("_HUMAN","",Protein)) %>%
               rename(prot=Protein)
    quadruple <- function(d,label) data.frame(Gene=label,chr=min(d$chr),start=min(d$start),end=max(d$end))
    caprion[c(55,237,433,435),]
    subset(glist_hg19,grepl("^AMY",gene))
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
    check_list <- c(55,237,278,385,390,433,435,845)
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
    write.csv(cis.vs.trans,file=file.path(analysis,"work",paste0("caprion",suffix,".cis.vs.trans")),row.names=FALSE,quote=FALSE)
    png(file.path(analysis,"work",paste0("caprion",suffix,".pqtl2d.png")),width=12,height=10,unit="in",res=300)
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
    htmlwidgets::saveWidget(r,file=file.path(analysis,"work",paste0("caprion",suffix,".pqtl2dplotly.html")))
    r <- qtl3dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                     snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                     gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                     TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                     xlab="pQTL position",ylab="Gene position")
    htmlwidgets::saveWidget(r,file=file.path(analysis,"work",paste0("caprion",suffix,".pqtl3dplotly.html")))
  '
}

function rsid()
# autosomes
{
  export cvt=${analysis}/work/caprion_dr.cis.vs.trans
  (
    echo snpid rsid
    sed '1d'  ${cvt} | cut -d, -f2 | sort | uniq | awk '{print "chr" $0}' | \
    join - ${INF}/work/INTERVAL.rsid | \
    sed 's/chr//'
  ) > ${analysis}/work/snpid_dr.lst
  (
    echo snpid rsid
    sed '1d'  ${cvt} | cut -d, -f2 | sort | uniq | \
    join -22 - <(cat ${analysis}/bgen/chr{1..22}.rsid | cut -d' ' -f3,4 | sort -k2,2)
  ) > ${analysis}/work/snpid_dr.lst2
}

for cmd in signals merge cistrans rsid; do $cmd; done
