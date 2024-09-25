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
    caprion <- dplyr::mutate(Protein, ProteinSNP=paste0(Protein,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
               dplyr::rename(prot=Protein,rsid=SNP,uniprot=Gene)
    caprion_dr <- dplyr::mutate(Protein_dr, ProteinSNP=paste0(Protein,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
                  dplyr::rename(prot=Protein,rsid=SNP,uniprot=Gene)
    b <- list()
    for(i in unique(dplyr::pull(caprion,chr)))
    {
       k <- dplyr::filter(caprion,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       dr <- dplyr::filter(caprion_dr,chr %in% i) %>%
             dplyr::select(chr,pos,uniprot,rsid,prot)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       b[[i]] <- pQTLtools::novelty_check(k,dr,ldops=list(bfile=bfile,plink=plink))
    }
    b[["23"]] <- dplyr::mutate(b[["X"]],known.seqnames="23",query.seqnames="23")
    r <- dplyr::filter(dplyr::bind_rows(b[-which(names(b)=="X")]),r2>=0.8) %>%
         dplyr::rename(known.gene=known.uniprot,query.gene=query.uniprot)
    s <- dplyr::mutate(r,seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(caprion,caprion_dr,r,s,file=file.path(analysis,"reports","caprion-caprion_dr.rda"))
    cvt <- "~/Caprion/analysis/work/caprion_dr.cis.vs.trans"
    f <- read.csv(cvt) %>%
         mutate(seqname=paste0("chr",SNPChrom),start=as.integer(SNPPos),end=SNPPos) %>%
         arrange(SNPChrom,SNPPos)
END
}

function deCODE()
{
  export v4=SomaLogicv4.tsv
  R --no-save -q <<END
    options(width=200)
    suppressMessages(library(dplyr))
    suppressMessages(require(openxlsx))
    analysis <- Sys.getenv("analysis")
    INF <- Sys.getenv("INF")
    plink <- Sys.getenv("plink")
    # Caprion
    load(file.path(analysis,"reports","caprion-caprion_dr.rda"))
    # deCODE
    xlsx <- file.path(Sys.getenv("deCODE"),"doc","ferkingstad21.xlsx")
    SomaLogicv4 <- openxlsx::read.xlsx(xlsx,sheet=1,startRow=3,colNames=TRUE,cols=1:12)
    out <- file.path(Sys.getenv("analysis"),"deCODE",Sys.getenv("v4"))
    write.table(SomaLogicv4,file=out,quote=FALSE,row.names=FALSE,sep="\t")
    sentinels <- openxlsx::read.xlsx(xlsx,sheet=2,startRow=3,colNames=TRUE) %>%
               select(10,9,3,7,11,14,15,18,19,27,29) %>%
               setNames(c("SeqId","uniprot","gene","prot","variant","chr","pos","A1","A2","beta","log10p")) %>%
               mutate(chr=gsub("chr","",chr),snpid=gap::chr_pos_a1_a2(chr,pos,A1,A2,prefix=""),
                      se=TwoSampleMR::get_se(beta,10^(-log10p)))
    rep38 <- sentinels %>%
             mutate(seqname=chr,start=as.integer(pos),end=pos) %>%
                    arrange(chr,pos)
    gr <- with(rep38,GenomicRanges::GRanges(seqnames=seqname,
                                    IRanges::IRanges(start,end,names=snpid),
                                    rsid=variant,chr38=chr,pos38=pos,prot=prot,uniprot=uniprot,
                                    beta=beta,se=se,log10p=log10p))
    suppressMessages(library(rtracklayer))
    path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- rtracklayer::import.chain(path)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    gr19 <- rtracklayer::liftOver(gr,ch) %>%
            data.frame %>%
            dplyr::mutate(seqnames=gsub("chr","",seqnames)) %>%
            dplyr::select(-group,-end,-width,-strand) %>%
            dplyr::rename(snpid=group_name,chr=seqnames,pos=start)
    write.table(gr19,file=file.path(analysis,"deCODE","deCODE.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
    suppressMessages(library(pQTLtools))
    b <- list()
    for(i in unique(dplyr::pull(caprion_dr,chr)))
    {
       dr <- dplyr::filter(caprion_dr,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       decode <- dplyr::filter(gr19,chr %in% i) %>%
                 dplyr::select(chr,pos,uniprot,snpid,prot) %>%
                 dplyr::rename(rsid=snpid)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       b[[i]] <- pQTLtools::novelty_check(decode,dr,ldops=list(bfile=bfile,plink=plink))
    }
    b[["23"]] <- dplyr::mutate(b[["X"]],known.seqnames="23",query.seqnames="23")
END
}

function UKB_PPP()
{
  R --no-save -q <<END
    options(width=2000)
    suppressMessages(require(dplyr))
    suppressMessages(require(openxlsx))
    suppressMessages(require(pQTLtools))
    analysis <- Sys.getenv("analysis")
    INF <- Sys.getenv("INF")
    plink <- Sys.getenv("plink")
    # Caprion
    load(file.path(analysis,"reports","caprion-caprion_dr.rda"))
    # UKB_PPP list
    results <- "/rds/project/jmmh2/rds-jmmh2-results/public/proteomics"
    url <- file.path(results,"UKB-PPP","doc","sun23.xlsx")
    ST10 <- read.xlsx(url,"ST10",startRow=3) %>%
            dplyr::mutate(uniprot=Target.UniProt,rsid=rsID,prot=Assay.Target) %>%
            dplyr::mutate(prot_rsid=paste0(uniprot,"-",rsid))
    caprion <- c(with(pQTLdata::caprion,uniprot),with(caprion_dr,uniprot)) %>%
               unique()
    overlap <- dplyr::filter(ST10,uniprot %in% caprion)
    dim(overlap)
    UKB_PPP <- dplyr::mutate(overlap,
               chrpos=strsplit(overlap[["Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)"]],":"),
               chr=as.integer(unlist(lapply(chrpos,"[[",1))),
               pos=as.integer(unlist(lapply(chrpos,"[[",2))),
               chrpos=paste(chr,pos,sep=":"))
    suppressMessages(require(GenomicRanges))
    replication <- read.delim(file.path(find.package("pQTLtools"),"tests","UKB-PPP.txt")) %>%
                   dplyr::select(known.seqnames,known.rsid,query.rsid,query.prot)
    variant_list <- unique(c(dplyr::pull(replication,known.rsid),
                             dplyr::pull(replication,query.rsid)))
    b <- list()
    for(i in unique(dplyr::pull(caprion_dr,Chromosome)))
    {
       u <- dplyr::filter(UKB_PPP,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       m <- dplyr::filter(caprion,Chromosome %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       b[[i]] <- novelty_check(u,m,ldops=list(bfile=bfile,plink=plink))
    }
    replication2 <- dplyr::filter(bind_rows(b), r2>=0.8)
    prot_rsid <- with(novel_data %>%
                 dplyr::left_join(gap.datasets::inf1[c("prot","gene")]),paste0(gene,"-",rsid))
    prot_rsid_repl <- with(replication2,paste0(query.prot,"-",query.rsid))
    novel <- setdiff(prot_rsid,prot_rsid_repl)
END
}

caprion
deCODE
UKB_PPP
