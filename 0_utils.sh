#!/usr/bin/bash

function pgwas()
{
  export protein=${1}
  export pqtl=${2}
# protein
  cd ~/Caprion/analysis/pgwas
  cat <(gunzip -c caprion-?-${protein}.fastGWA.gz | head -1 | paste <(echo Batch-Protein) -|sed 's/-/\t/') \
      <(zgrep -w ${pqtl} caprion-?-${protein}.fastGWA.gz) | \
  sed 's/caprion-//;s/.fastGWA.gz:/\t/;s/-/\t/' | \
  Rscript -e '
    suppressMessages(library(dplyr))
    pqtl <- Sys.getenv("pqtl")
    d <- read.table("stdin",header=TRUE) %>%
         arrange(Batch)
    knitr::kable(d,caption=paste("Effect sizes of",pqtl),digits=3)
  '
  cat <(gunzip -c ~/Caprion/analysis/METAL/${protein}-1.tbl.gz | head -1 | paste <(echo Protein) -) \
      <(zgrep -H -w ${pqtl} ~/Caprion/analysis/METAL/${protein}-1.tbl.gz) | \
  sed 's/INHBE//;s/-1.tbl.gz:/\t/;s/:/\t/' | \
  Rscript -e '
    suppressMessages(library(dplyr))
    pqtl <- Sys.getenv("pqtl")
    d <- read.table("stdin",header=TRUE) %>%
         arrange(Protein)
    knitr::kable(d,caption=paste("Effect sizes of",pqtl),digits=3)
  '
  cat <(gunzip -c ~/Caprion/analysis/METAL3/${protein}-1.tbl.gz | head -1 | paste <(echo Protein) -) \
      <(zgrep -H -w ${pqtl} ~/Caprion/analysis/METAL3/${protein}-1.tbl.gz) | \
  sed 's/INHBE//;s/-1.tbl.gz:/\t/;s/:/\t/' | \
  Rscript -e '
    suppressMessages(library(dplyr))
    pqtl <- Sys.getenv("pqtl")
    d <- read.table("stdin",header=TRUE) %>%
         arrange(Protein)
    knitr::kable(d,caption=paste("Effect sizes of",pqtl),digits=3)
  '
# peptides
  cd ~/Caprion/analysis/peptide/${protein}/
  cat <(gunzip -c ${protein}-?-*.fastGWA.gz | head -1 | paste <(echo Batch-Peptide) -|sed 's/-/\t/') \
      <(zgrep -w ${pqtl} *fastGWA.gz) | \
  sed 's/INHBE-//;s/.fastGWA.gz:/\t/;s/-/\t/' | \
  Rscript -e '
    suppressMessages(library(dplyr))
    pqtl <- Sys.getenv("pqtl")
    d <- read.table("stdin",header=TRUE) %>%
         arrange(Peptide,Batch)
    knitr::kable(d,caption=paste("Effect sizes of",pqtl),digits=3)
  '
  cat <(gunzip -c METAL/*-1.tbl.gz | head -1 | paste <(echo Peptide) -) \
      <(zgrep -H -w ${pqtl} METAL/*-1.tbl.gz) | \
  sed 's|METAL/||;s/-1.tbl.gz:/\t/;s/:/\t/' | \
  Rscript -e '
    suppressMessages(library(dplyr))
    pqtl <- Sys.getenv("pqtl")
    d <- read.table("stdin",header=TRUE) %>%
         arrange(Peptide)
    knitr::kable(d,caption=paste("Effect sizes of",pqtl),digits=3)
  '
}


function fp_data()
{
  cp  ${analysis}/work/caprion${suffix}.merge ${analysis}/work/tbl.tsv
  cut -f2-4 ${analysis}/work/tbl.tsv | \
  awk 'NR>1' | \
  sort -k1,2n | \
  uniq | \
  awk -vOFS="\t" '{print $1":"$2,$3}' > ${analysis}/work/rsid.tsv
  (
    gunzip -c ${analysis}/pgwas${suffix}/caprion-*fastGWA.gz | head -1
    awk 'NR>1' ${analysis}/work/tbl.tsv | \
    cut -f1,4,14 --output-delimiter=' ' | \
    parallel -j10 -C' ' '
      export direction=$(zgrep -w {2} ${analysis}/METAL${suffix}/{1}${suffix}-1.tbl.gz | cut -f13)
      let j=1
      for i in $(grep "Input File" ${analysis}/METAL${suffix}/{1}${suffix}-1.tbl.info | cut -d" " -f7)
      do
         export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
         if [ "$n" != "?" ]; then zgrep -H -w {2} $i; fi
         let j=$j+1
      done
  '
  ) | \
  sed 's/.gz//g' > ${analysis}/work/all.tsv
}

function fp()
{
  if [ ! -d ${analysis}/METAL${suffix}/fp ]; then mkdir -o mkdir -p ${analysis}/METAL${suffix}/fp; fi
  Rscript -e '
    require(gap)
    require(dplyr)
    suffix <- Sys.getenv("suffix")
    tbl <- read.delim("~/Caprion/analysis/work/tbl.tsv") %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position)) %>%
           arrange(prot,SNP)
    all <- read.delim("~/Caprion/analysis/work/all.tsv") %>%
           rename(EFFECT_ALLELE=A1,REFERENCE_ALLELE=A2) %>%
           mutate(CHR=gsub(suffix,"",CHR),
                  CHR=gsub("/home/jhz22/Caprion/analysis/pgwas/caprion-|.fastGWA","",CHR)) %>%
           mutate(batch_prot_chr=strsplit(CHR,"-|:"),
                  batch=unlist(lapply(batch_prot_chr,"[",1)),
                  prot=unlist(lapply(batch_prot_chr,"[",2)),
                  CHR=unlist(lapply(batch_prot_chr,"[",3)),CHR=gsub("chrX","23",CHR)) %>%
           mutate(MarkerName=paste0(CHR,":",POS),
                  study=case_when(batch == "1" ~ paste0("1. ZWK (",N,")"),
                                  batch == "2" ~ paste0("2. ZYQ (",N,")"),
                                  batch == "3" ~ paste0("3. UDP (",N,")"),
                                  TRUE ~ "---")) %>%
           arrange(study) %>%
           select(-batch_prot_chr)
    rsid <- read.table("~/Caprion/analysis/work/rsid.tsv",col.names=c("MarkerName","rsid"))
    pdf("~/Caprion/analysis/work/fp.pdf",width=8,height=5)
    METAL_forestplot(tbl,all,rsid,package="metafor",method="FE",cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="Effect")
    dev.off()
  '
}

function HetISq()
# Code extracted from caprion.Rmd
{
  Rscript -e '
    suppressMessages(require(dplyr))
    suffix <- Sys.getenv("suffix")
    all <- read.delim("~/Caprion/analysis/work/all.tsv") %>%
           mutate(CHR=gsub("/home/jhz22/Caprion/analysis/work/pgwas/caprion-|.fastGWA","",CHR),CHR=gsub(suffix,"",CHR)) %>%
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
    b <- full_join(b1,b2,by=c('prot'='prot.ZYQ','SNP'='SNP.ZYQ')) %>%
         full_join(b3,by=c('prot'='prot.UDP','SNP'='SNP.UDP')) %>%
         mutate(directions=gsub("NA","?",paste0(direction.ZWK,direction.ZYQ,direction.UDP))) %>%
         select(-Batch.ZWK,-Batch.ZYQ,-Batch.UDP,direction.ZWK,direction.ZYQ,direction.UDP)
    tbl <- read.delim("~/Caprion/analysis/work/tbl.tsv") %>%
           arrange(prot,MarkerName) %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position), index=1:n())
    Het <- filter(tbl,HetISq>=75) %>%
           select(prot,SNP,Direction,HetISq,index) %>%
           left_join(select(b,prot,SNP,P.ZWK,P.ZYQ,P.UDP,BETA.ZWK,BETA.ZYQ,BETA.UDP))
    write.csv(Het,file=file.path(analysis,"work", paste0("HetISq75",suffix,".csv"),row.names=FALSE,quote=FALSE)
    write(Het[['index']],file=file.path(analysis, "work", paste0("HetISq75",suffix,".index"),sep=",",ncolumns=nrow(Het))
  '
}

function fplz()
{
  export metal=${analysis}/METAL${suffix}
# HSPB1_rs114800762 is missing as dug by the following code.
  join -a1 <(sed '1d' ${analysis}/work/caprion${suffix}.merge | awk '{print $1"_"$4}' | sort -k1,1 ) \
           <(ls ${analysis}/METAL${suffix}/qqmanhattanlz/lz/*pdf | xargs -l basename -s .pdf | awk '{print $1,NR}') | \
  awk 'NF<2' | \
  sed 's/_/ /' | \
  parallel -C' ' 'ls ${analysis}/METAL${suffix}/qqmanhattanlz/lz/{1}*pdf'
# forest/locuszoom left-right format
  ulimit -n
  ulimit -S -n 2048
  qpdf --empty --pages $(sed '1d' ${analysis}/work/caprion${suffix}.merge | sort -k1,1 -k4,4 | cut -f1,4 --output-delimiter=' ' | \
                         parallel -C' ' 'ls $(echo ${analysis}/METAL${suffix}/qqmanhattanlz/lz/{1}_{2}.pdf | sed "s/:/_/")') -- lz2.pdf
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
  qpdf fp+lz.pdf --pages . \
                 $(sed '1d' ${analysis}/work/caprion${suffix}.merge | sort -k1,1 -k4,4 | awk '$15>=75{printf " "NR}' | sed 's/ //;s/ /,/g') \
       -- HetISq75.pdf
}

function pdf()
{
  export f=${analysis}/work/caprion${suffix}.signals
  export N=$(sed '1d' ${f} | wc -l)
  export g=10
  export d=${analysis}/METAL${suffix}/qqmanhattanlz
  module load ceuadmin/pdfjam gcc/6
# qq-manhattan
  ls *_qq.png | xargs -l basename -s _qq.png | \
  parallel -C' ' 'convert -resize 150% {}_qq.png {}_qq.pdf;convert {}_manhattan.png {}_manhattan.pdf'
  qpdf --empty --pages $(ls *_qq.pdf) -- qq.pdf
  qpdf --empty --pages $(ls *_manhattan.pdf) -- manhattan.pdf
  pdfseparate qq.pdf temp-%04d-qq.pdf
  pdfseparate manhattan.pdf temp-%04d-manhattan.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR<=500') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile qq-manhattan1.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR>500 && NR<=1000') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile qq-manhattan2.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR>1000 && NR<=1500') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile qq-manhattan3.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR>1500') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile qq-manhattan4.pdf
  qpdf --empty --pages qq-manhattan*pdf -- qq-manhattan.pdf
  rm temp*
# lz
  sed '1d' ${f} | \
  awk -vN=${N} -vg=${g} '
  function ceil(v) {return(v+=v<0?0:0.999)}
  {
     gsub(":","_",$7)
     printf "%d %s %d %d %s\n", ceil(NR*g/N), $1, $2, $4, $7
  } ' > ${N}
  for i in `seq ${g}`
  do
     export n=$(awk -v i=${i} '$1==i' ${N} | wc -l)
     export n2=$(expr ${n} \* 2)
     qpdf --empty --pages $(awk -v i=${i} '$1==i' ${N} | \
                            awk -v d=${d} -v suffix=${suffix} '{print d"/"$2 suffix"_"$5".pdf"}' | \
                            sort -k1,1 | \
                            tr '\n' ' ';echo) \
          -- lz2-${i}.pdf
     qpdf --pages . 1-${n2}:odd -- lz2-${i}.pdf lz-${i}.pdf
     rm lz2-${i}.pdf
  done
  qpdf --empty --pages $(echo lz-{1..10}.pdf) -- lz.pdf
  rm ${N}
# fp-lz
  pdfseparate ${analysis}/work/fp.pdf temp-%04d-fp.pdf
  pdfseparate ${analysis}/METAL${suffix}/qqmanhattanlz/lz.pdf temp-%04d-lz.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR<=500') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile fp-lz1.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR>500 && NR<=1000') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile fp-lz2.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR>1000 && NR<=1500') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile fp-lz3.pdf
  pdfjam $(ls temp-*-*.pdf|awk 'NR>1500') --nup 2x1 --landscape --papersize '{5in,16in}' --outfile fp-lz4.pdf
  qpdf --empty --pages fp-lz*pdf -- fp-lz.pdf
  rm temp*
}

function INHBE
{
(
  pgwas INHBE rs149830883
  pgwas INHBE rs11172187
) > ~/Caprion/analysis/work/INHBE.txt
}

function hist_corr_lm()
{
Rscript -e '
  options(width=200)
  suppressMessages(library(Biobase))
  suppressMessages(library(dplyr))
  suppressMessages(library(Hmisc))
  suppressMessages(library(pheatmap))
  filter(pQTLdata::caprion,Protein=="INHBE_HUMAN") %>%
  select(Protein,Accession,Gene,Protein.Description)
  protein_peptide <- function(protein="INHBE",suffix="ZWK")
  {
    cat("\n**",suffix,"**\n",sep="")
    dir <- "~/rds/projects/Caprion_proteomics"
    load(file.path(dir,"pilot",paste(suffix,"rda",sep=".")))
    n <- paste("protein",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% paste(protein,"HUMAN",sep="_")
    prot <- matrix(p[g,],nrow=1,dimnames=list(protein,names(p[g,]))) %>%
            data.frame
    n <- paste("dr",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% paste(protein,"HUMAN",sep="_")
    dr <- matrix(p[g,],nrow=1,dimnames=list(protein,names(p[g,]))) %>%
          data.frame
    n <- paste("mapping",suffix,sep="_")
    m <- subset(get(n),grepl(protein,Protein))
    igID <- m[["Isotope.Group.ID"]]
    n <- paste("peptide",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% igID
    pept <- data.frame(p[g,])
    prot_pept <- bind_rows(pept,prot)
    n <- rownames(prot_pept)
    d <- t(prot_pept) %>%
         data.frame() %>%
         setNames(n)
    cat("\nProtein/Protein_DR")
    s1 <- summary(lm(INHBE~d[["442593377"]]+d[["442626845"]]+d[["442628596"]]+d[["442658425"]],data=d))
    dr_pept <- bind_rows(pept,dr)
    n <- rownames(dr_pept)
    d <- t(dr_pept) %>%
         data.frame() %>%
         setNames(n)
    s2 <- summary(lm(INHBE~d[["442593377"]]+d[["442626845"]]+d[["442628596"]]+d[["442658425"]],data=d))
    print(knitr::kable(cbind(coef(s1),coef(s2)),digits=3))
    opar <- par()
    png(file.path(dir,"analysis","work",paste(protein,suffix,"dist.png",sep="-")),
        width=12,height=10,units="in",pointsize=4,res=300)
    par(mar=c(15,10,5,5), font=2, font.lab = 5, font.axis = 5)
    source("https://raw.githubusercontent.com/jinghuazhao/tests/main/Hmisc/hist.data.frame.R")
    hist.data.frame(d,cex.axis=2.5,cex.lab=2.5,cex.mtext=2.5,cex.names=2.5,ylab=expression("Frequency"))
    dev.off()
    par(opar)
    png(file.path(dir,"analysis","work",paste(protein,suffix,"corr.png",sep="-")),
        width=12,height=10,units="in",pointsize=4,res=300)
    pheatmap(cor(t(prot_pept)),display_numbers=TRUE,fontsize=24)
    dev.off()
    write.csv(m[c("Isotope.Group.ID","Modified.Peptide.Sequence","Protein")],
              file=file.path(dir,"analysis","work",paste(protein,suffix,"mapping.csv",sep="-")),
              quote=FALSE,row.names=FALSE)
    prot_pept
  }
  zwk <- protein_peptide()
  zyq <- protein_peptide(suffix="ZYQ")
  udp <- protein_peptide(suffix="UDP")
'
}

function barplot()
{
Rscript -e '
  one <- read.delim("1")
  batches <- unique(with(one,Batch))
  peptides <- unique(with(one,Peptide))

  m <- s <- matrix(NA,length(peptides),length(batches))
  colnames(m) <- paste(batches)
  colnames(s) <- paste(batches)
  rownames(m) <- paste(peptides)
  rownames(s) <- paste(peptides)
  for(p in paste(peptides)) for(b in paste(batches))
  {
    d <- subset(one, Peptide==p & Batch==b)
    m[p,b] <- d[["BETA"]]
    s[p,b] <- d[["SE"]]
  }

  s.bar <- function(x, y, upper, lower=upper, length=0.1,...)
  {
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }

  png("INHBE-peptides.png",res=300,width=6,height=6,units="in")
  z <- barplot(m , beside=TRUE , legend.text=TRUE, args.legend=c(x=6,y=-0.9),
               col=c("blue" , "skyblue", "red", "green"), xlab="Batch", ylab="Beta", ylim=c(-1.2,0))
  s.bar(z,m,s)
  title("INHBE rs149830883 association")
  dev.off()
'
}

export analysis=~/Caprion/analysis
export suffix=_dr

function ucsc_annotate()
{
  Rscript -e '
    options(width=200)
    library(dplyr)
    analysis <- Sys.getenv("analysis")
    suffix <- Sys.getenv("suffix")
    library(pQTLdata)
    nodup <- function(x) sapply(x, function(s) unique(unlist(strsplit(s,";")))[1])
    ucsc <- hg19Tables %>%
            group_by(acc) %>%
            summarize(
                 prot=paste(uniprotName,collapse=";"),
                 chrom=paste(X.chrom,collapse=";"),
                 start=min(chromStart),
                 end=max(chromEnd),
                 gene=paste(geneName,collapse=";")
            )
    # uniprot IDs are the same if proteins are the same
    p <- select(caprion,Accession,Protein,Gene) %>%
         left_join(ucsc,by=c("Protein"="prot")) %>%
         select(Accession,Protein,Gene,gene,acc,chrom,start,end)
    # however even with same uniprotID their protein names may be different
    u <- select(caprion,Accession,Protein,Gene) %>%
         left_join(ucsc,by=c("Accession"="acc")) %>%
         mutate(chrom=nodup(chrom)) %>%
         filter(!is.na(Protein)) %>%
         select(Accession,Protein,gene,Gene,prot,chrom,start,end)
    # The following check shows merge by uniprot is more sensible
    filter(p,Accession!=acc)
    filter(p,Gene!=gene)
    filter(p,is.na(start))
    filter(u,Protein!=prot)
    umiss <- with(u,is.na(start))
    filter(u,umiss) %>% pull(Accession)
    # "P04745" "P02655" "P55056" "P0C0L5" "P62805" "P69905"
    # confirmed form UniProt.org that Gene is more up-to-date
    # (obsolute), (APOC2), (APOC4), (C4B; C4B_2), (H4C1; H4C2; H4C3; H4C4; H4C5; H4C6; H4C8; H4C9; H4C11; H4C12; H4C13; H4C14; H4C15; H4C16), (HBA1; HBA2)
    # They are amended according to glist-hg19 in the function below.
    u[umiss,"Protein"] <- paste0(c("AMY1","APOC2","APOC4","CO4B","H4","HBA"),"_HUMAN")
    u[umiss,"Gene"] <- c("AMY1","APOC2","APOC4","CO4B","H4C","HBA")
    u[umiss,"chrom"] <- c("chr1","chr19","chr19","chr6","chr6","chr16")
    u[umiss,"start"] <- c(104198140,45449238,45445494,31949833,26021906,222845)
    u[umiss,"end"] <- c(104301311,45452822,45452822,32003195,27841289,227520)
    caprion_modified <- u
    a <- filter(u,umiss) %>%
         transmute(acc=Accession,prot=Protein,gene=Gene,chrom,start,end)
    ucsc2 <- ucsc %>%
             mutate(prot=nodup(prot),chrom=nodup(chrom),gene=nodup(gene)) %>%
             mutate(chrom=gsub("chrX","chr23",chrom),chrom=gsub("chrY","chr24",chrom))
             bind_rows(a)
    load("~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda")
    refgene_gene_coordinates_h19 <- ucsc2 %>%
                                    transmute(chromosome=gsub("chr","",chrom),
                                              gene_transcription_start=start,
                                              gene_transcription_stop=end,
                                              gene_name=gene,acc,prot,
                                              gene_transcription_midposition=(start+end)/2)
    save(ld_block_breaks_pickrell_hg19_eur,refgene_gene_coordinates_h19,file="ucsc_hg19_reference_data.rda")
    cis.vs.trans <- read.csv(file=file.path(analysis,"work",paste0("caprion",suffix,".cis.vs.trans"))) %>%
                    arrange(prot,SNPChrom,SNPPos,Type) %>%
                    transmute(prot,chrom=paste0("chr",SNPChrom),start=SNPPos,end=SNPPos,cistrans=Type)
    library(valr)
    d <- bed_intersect(cis.vs.trans,ucsc2) %>%
         transmute(chromosome=gsub("chr","",chrom),position=start.x,nearest_gene_name=gene.y,cistrans=cistrans.x,protein=prot.x)
    write.table(d,file="~/cambridge-ceu/turboman/caprion.txt",quote=FALSE,row.names=FALSE)
'
}

function vep_annotate()
{
  if [ ! -d ${analysis}/METAL${suffix}/vep ]; then mkdir ${analysis}/METAL${suffix}/vep; fi
  export cvt=${analysis}/work/caprion${suffix}.cis.vs.trans
  cut -d"," -f3 ${cvt} | \
  sort -k1,1 | \
  uniq | \
  parallel -C' ' '
    (
      echo "##fileformat=VCFv4.0"
      echo "#CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO"
      awk -vFS="," -vprot={} "\$3==prot {print \$2}" ${cvt} | \
      sort -k1,1 | \
      zgrep -f - -w ${analysis}/METAL${suffix}/{}${suffix}-1.tbl.gz | \
      cut -f1-5 | \
      awk "{print \$1,\$2,\$3,toupper(\$4),toupper(\$5),\".\",\".\",\".\"}"
    ) | \
    tr " " "\t" > ${analysis}/METAL${suffix}/vep/{}${suffix}.vcf
  # VEP annotation
    cd ${HPC_WORK}/loftee
    vep --input_file ${analysis}/METAL${suffix}/vep/{1}${suffix}.vcf \
        --output_file ${analysis}/METAL${suffix}/vep/{1}${suffix}.tab --force_overwrite \
        --cache --dir_cache ${HPC_WORK}/ensembl-vep/.vep --dir_plugins ${HPC_WORK}/loftee --offline \
        --species homo_sapiens --assembly GRCh37 --pick --nearest symbol --symbol --plugin TSSDistance \
        --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz \
        --tab
    cd -
    (
      echo chromosome position nearest_gene_name cistrans
      awk -vFS="," -vprot={} "\$3==prot {print \$2,\$8,\$9,\$10}" ${cvt} | \
      sort -k1,1 | \
      join - <(awk "!/#/{print \$1,\$21}" ${analysis}/METAL${suffix}/vep/{}${suffix}.tab | sort -k1,1) | \
      awk "{print \$2,\$3,\$5,\$4}" | \
      sort -k1,1n -k2,2n
    ) > ${analysis}/METAL${suffix}/vep/{}${suffix}.txt
  '
}

function signal_comparison()
{
  Rscript -e '
    options(width=200)
    library(dplyr)
    cvt <- read.csv("~/Caprion/analysis/work/caprion.cis.vs.trans") %>%
           arrange(prot) %>%
           mutate(chrom=paste0("chr",SNPChrom),start=SNPPos,end=SNPPos)
    dim(cvt)
    head(cvt)
    cvt_dr <- read.csv("~/Caprion/analysis/work/caprion_dr.cis.vs.trans") %>%
              arrange(prot) %>%
              mutate(chrom=paste0("chr",SNPChrom),start=SNPPos,end=SNPPos)
    dim(cvt_dr)
    head(cvt_dr)
    library(valr)
    intersect(cvt,cvt_dr)
    right_join(cvt,cvt_dr,by=c("prot","SNP"))
    intersect(select(cvt,chrom,start,end),select(cvt_dr,chrom,start,end)) %>% nrow
  # 394
    intersect(select(cvt,chrom,start,end,prot,SNP),select(cvt_dr,chrom,start,end,prot,SNP)) %>% dim
  # 446
  # potential to add novelty check
  '
}

function ukb_ppp()
{
  export rt=~/rds/results/public/proteomics/UKB-PPP/sun23
  export f=A1BG.tsv
  gunzip -c ${rt}/UKB-PPP\ pGWAS\ summary\ statistics\ \(reformatted\)/European\ \(discovery\)/A1BG_*gz | \
  awk 'NR==1||$13>=7.30103' > ${f}
  Rscript -e '
    options(width=200)
    library(dplyr)
    f <- Sys.getenv("f")
    d <- read.delim(f)
    d  <- within(d,{LOG10P=-LOG10P})
    qtls <- qtlFinder(d,Chromosome="CHROM",Position="GENPOS",
                      MarkerName="ID",Allele1="ALLELE0",Allele2="ALLELE1",
                      EAF="A1FREQ",Effect="BETA",StdErr="SE",log10P="LOG10P") %>%
            mutate(rsid=gsub(":imp:v1","",rsid)) %>%
            select(-a1,-a2,-.overlap)
  '
}
