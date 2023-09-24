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
    tr " " "\t" > ${analysis}/METAL${suffix}/vep/{}.vcf
  # VEP annotation
    cd ${HPC_WORK}/loftee
    vep --input_file ${analysis}/METAL${suffix}/vep/{1}.vcf --output_file ${analysis}/METAL${suffix}/vep/{1}.tab --force_overwrite \
        --cache --dir_cache ${HPC_WORK}/ensembl-vep/.vep --dir_plugins ${HPC_WORK}/loftee --offline \
        --species homo_sapiens --assembly GRCh37 --pick --nearest symbol --symbol --plugin TSSDistance \
        --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz \
        --tab
    cd -
    awk -vFS="," -vprot={} "\$3==prot {print \$2,\$8,\$9,\$10}" ${cvt} | \
    sort -k1,1 | \
    join - <(awk "!/#/{print \$1,\$21}" ${analysis}/METAL${suffix}/vep/{}.tab | sort -k1,1) | \
    awk "{print \$2,\$3,\$5,\$4}" | \
    sort -k1,1n -k2,2n > ${analysis}/METAL${suffix}/vep/{}.txt
  '
}

vep_annotate
