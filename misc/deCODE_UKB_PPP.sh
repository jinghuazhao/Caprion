#!/usr/bin/bash

export deCODE=~/rds/results/public/proteomics/deCODE
export analysis=~/Caprion/analysis
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export plink=/rds/user/jhz22/hpc-work/bin/plink

if [ ! -d ${analysis}/deCODE ]; then mkdir -p ${analysis}/deCODE; fi
if [ ! -d ${analysis}/UKB_PPP ]; then mkdir -p ${analysis}/UKB_PPP; fi

function caprion()
{
  R --no-save -q <<END
    suppressMessages(library(dplyr))
    suppressMessages(library(gap))
    require(stringr)
    options(width=2000)
    analysis <- Sys.getenv("analysis")
    INF <- Sys.getenv("INF")
    plink <- Sys.getenv("plink")
    freq <- read.table("~/Caprion/analysis/bgen/caprion.freq",col.names=c("SNP","REF","ALT","ALT_FREQS"))
    tbl <- read.delim(file.path(analysis,"work","tbl.tsv")) %>%
           dplyr::select(prot,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,Direction,HetISq,logHetP,N) %>%
           dplyr::mutate(Allele1=toupper(Allele1),Allele2=toupper(Allele2)) %>%
           dplyr::rename(Protein=prot,SNP=MarkerName,EA=Allele1,OA=Allele2,EAF=Freq1,Log10P=log.P.)
    tbl_dr <- read.delim(file.path(analysis,"work","tbl_dr.tsv")) %>%
              dplyr::select(prot,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,Direction,HetISq,logHetP,N) %>%
              dplyr::mutate(Allele1=toupper(Allele1),Allele2=toupper(Allele2)) %>%
              dplyr::rename(Protein=prot,SNP=MarkerName,EA=Allele1,OA=Allele2,EAF=Freq1,Log10P=log.P.)
    Protein <- read.csv(file.path(analysis,"work","caprion.cis.vs.trans")) %>%
               dplyr::mutate(geneRegion=paste(geneChrom,":",geneStart,"-",geneEnd)) %>%
               dplyr::left_join(freq) %>%
               dplyr::arrange(prot,SNPChrom,SNPPos) %>%
               dplyr::rename(Protein=prot) %>%
               dplyr::select(-log10p,-geneChrom,-geneStart,-geneEnd,-SNPChrom,-SNPPos,-cis) %>%
               dplyr::left_join(tbl) %>%
               dplyr::select(Gene,Protein,geneRegion,SNP,Type,EA,OA,EAF,Effect,StdErr,Log10P,Direction,HetISq,logHetP,N,REF,ALT,ALT_FREQS)
    Protein_dr <- read.csv(file.path(analysis,"work","caprion_dr.cis.vs.trans")) %>%
               dplyr::mutate(geneRegion=paste(geneChrom,":",geneStart,"-",geneEnd)) %>%
               dplyr::left_join(freq) %>%
               dplyr::arrange(prot,SNPChrom,SNPPos) %>%
               dplyr::rename(Protein=prot) %>%
               dplyr::select(-log10p,-geneChrom,-geneStart,-geneEnd,-SNPChrom,-SNPPos,-cis) %>%
               dplyr::left_join(tbl_dr) %>%
               dplyr::select(Gene,Protein,geneRegion,SNP,Type,EA,OA,EAF,Effect,StdErr,Log10P,Direction,HetISq,logHetP,N,REF,ALT,ALT_FREQS)
    panel <- pQTLdata::caprion %>%
             transmute(prot=gsub("_HUMAN","",Protein),uniprot=Accession)
    caprion <- dplyr::mutate(Protein, ProteinSNP=paste0(Protein,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
               dplyr::rename(prot=Protein,rsid=SNP) %>%
               dplyr::left_join(panel)
    caprion_dr <- dplyr::mutate(Protein_dr, ProteinSNP=paste0(Protein,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
                  dplyr::rename(prot=Protein,rsid=SNP) %>%
                  dplyr::left_join(panel)
    z <- list()
    for(i in unique(dplyr::pull(caprion,chr)))
    {
       k <- dplyr::filter(caprion,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       dr <- dplyr::filter(caprion_dr,chr %in% i) %>%
             dplyr::select(chr,pos,uniprot,rsid,prot)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       z[[i]] <- pQTLtools::novelty_check(k,dr,ldops=list(bfile=bfile,plink=plink))
    }
    z[["23"]] <- dplyr::mutate(z[["X"]],known.seqnames="23",query.seqnames="23")
    r <- dplyr::filter(dplyr::bind_rows(z[-which(names(z)=="X")]),r2>=0.8) %>%
         dplyr::rename(known.gene=known.uniprot,query.gene=query.uniprot)
    s <- dplyr::mutate(r,seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(panel,caprion,caprion_dr,r,s,file=file.path(analysis,"reports","caprion-caprion_dr.rda"))
    cvt <- "~/Caprion/analysis/work/caprion_dr.cis.vs.trans"
    f <- read.csv(cvt) %>%
         mutate(seqname=paste0("chr",SNPChrom),start=as.integer(SNPPos),end=SNPPos) %>%
         arrange(SNPChrom,SNPPos)
END
}

function deCODE()
{
  R --no-save -q <<END
    options(width=200)
    suppressMessages(library(dplyr))
    suppressMessages(require(openxlsx))
    suppressMessages(library(rtracklayer))
    suppressMessages(library(pQTLtools))
    suppressMessages(require(GenomicRanges))
    analysis <- Sys.getenv("analysis")
    INF <- Sys.getenv("INF")
    plink <- Sys.getenv("plink")
    # Caprion
    load(file.path(analysis,"reports","caprion-caprion_dr.rda"))
    # deCODE
    xlsx <- file.path(Sys.getenv("deCODE"),"doc","ferkingstad21.xlsx")
    SomaLogicv4 <- openxlsx::read.xlsx(xlsx,sheet=1,startRow=3,colNames=TRUE,cols=1:12)
    sentinels <- openxlsx::read.xlsx(xlsx,sheet=2,startRow=3,colNames=TRUE) %>%
               dplyr::select(10,9,3,7,11,14,15,18,19,27,29) %>%
               setNames(c("SeqId","uniprot","gene","prot","variant","chr","pos","A1","A2","beta","log10p")) %>%
               dplyr::mutate(chr=gsub("chr","",chr),snpid=gap::chr_pos_a1_a2(chr,pos,A1,A2,prefix=""),
                             se=TwoSampleMR::get_se(beta,10^(-log10p)))
    rep38 <- sentinels %>%
             dplyr::mutate(seqname=chr,start=as.integer(pos),end=pos) %>%
             dplyr::arrange(chr,pos)
    gr <- with(rep38,GenomicRanges::GRanges(seqnames=seqname,
                                    IRanges::IRanges(start,end,names=snpid),
                                    rsid=variant,chr38=chr,pos38=pos,prot=prot,uniprot=uniprot,
                                    beta=beta,se=se,log10p=log10p))
    f <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- rtracklayer::import.chain(f)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    deCODE <- rtracklayer::liftOver(gr,ch) %>%
              data.frame %>%
              dplyr::mutate(seqnames=gsub("chr","",seqnames)) %>%
              dplyr::select(-group,-end,-width,-strand) %>%
              dplyr::rename(snpid=group_name,chr=seqnames,pos=start)
    z <- list()
    for(i in unique(dplyr::pull(caprion_dr,chr)))
    {
       dr <- dplyr::filter(caprion_dr,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       decode <- dplyr::filter(deCODE,chr %in% i) %>%
                 dplyr::select(chr,pos,uniprot,snpid,prot) %>%
                 dplyr::rename(rsid=snpid)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       z[[i]] <- pQTLtools::novelty_check(decode,dr,ldops=list(bfile=bfile,plink=plink))
    }
    z[["23"]] <- dplyr::mutate(z[["X"]],known.seqnames="23",query.seqnames="23")
    u <- dplyr::filter(dplyr::bind_rows(z[-which(names(z)=="X")]),r2>=0.8) %>%
         dplyr::rename(known.gene=known.uniprot,query.gene=query.uniprot)
    v <- dplyr::mutate(u,seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(deCODE,u,v,file=file.path(analysis,"deCODE","deCODE.rda"))
    f <- file.path(Sys.getenv("analysis"),"deCODE","SomaLogicv4.tsv")
    write.table(SomaLogicv4,file=f,quote=FALSE,row.names=FALSE,sep="\t")
    write.table(gr19,file=file.path(analysis,"deCODE","deCODE.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
END
}

function UKB_PPP()
{
  R --no-save -q <<END
    options(width=2000)
    suppressMessages(require(dplyr))
    suppressMessages(require(openxlsx))
    suppressMessages(require(pQTLtools))
    suppressMessages(require(GenomicRanges))
    analysis <- Sys.getenv("analysis")
    INF <- Sys.getenv("INF")
    plink <- Sys.getenv("plink")
    # Caprion
    load(file.path(analysis,"reports","caprion-caprion_dr.rda"))
    # UKB_PPP list
    results <- "/rds/project/jmmh2/rds-jmmh2-results/public/proteomics"
    url <- file.path(results,"UKB-PPP","doc","sun23.xlsx")
    ST10 <- read.xlsx(url,"ST10",startRow=4) %>%
            dplyr::mutate(uniprot=Target.UniProt,rsid=rsID,prot=Assay.Target) %>%
            dplyr::mutate(prot_rsid=paste0(uniprot,"-",rsid))
    uniprots <- c(with(panel,uniprot),with(caprion_dr,uniprot)) %>%
                unique()
    overlap <- dplyr::filter(ST10,uniprot %in% uniprots)
    dim(overlap)
    UKB_PPP <- dplyr::mutate(overlap,
               chrpos=strsplit(overlap[["Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)"]],":"),
               chr=unlist(lapply(chrpos,"[[",1)),
               pos=unlist(lapply(chrpos,"[[",2)),
               a1=unlist(lapply(chrpos,"[[",3)),
               a2=unlist(lapply(chrpos,"[[",4))) %>%
               dplyr::filter(!is.na(a1)&!is.na(a2)) %>%
               dplyr::mutate(rsid=gap::chr_pos_a1_a2(chr,pos,a1,a2)) %>%
               dplyr::select(-chrpos)
    z <- list()
    for(i in unique(dplyr::pull(caprion_dr,chr)))
    {
       u <- dplyr::filter(UKB_PPP,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       m <- dplyr::filter(caprion_dr,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       z[[i]] <- novelty_check(u,m,ldops=list(bfile=bfile,plink=plink))
    }
    z[["23"]] <- dplyr::mutate(z[["X"]],known.seqnames="23",query.seqnames="23")
    w <- dplyr::filter(dplyr::bind_rows(z[-which(names(z)=="X")]),r2>=0.8) %>%
         dplyr::rename(known.gene=known.uniprot,query.gene=query.uniprot)
    x <- dplyr::mutate(w,seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(UKB_PPP,w,x,file=file.path(analysis,"deCODE","deCODE.rda"))
    prot_rsid <- with(novel_data %>%
                 dplyr::left_join(gap.datasets::inf1[c("prot","gene")]),paste0(gene,"-",rsid))
    prot_rsid_repl <- with(replication2,paste0(query.prot,"-",query.rsid))
    novel <- setdiff(prot_rsid,prot_rsid_repl)
END
}

caprion
deCODE
UKB_PPP
