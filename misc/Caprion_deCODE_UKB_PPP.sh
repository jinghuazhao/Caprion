#!/usr/bin/bash

export deCODE=~/rds/results/public/proteomics/deCODE
export analysis=~/Caprion/analysis
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export plink=/rds/user/jhz22/hpc-work/bin/plink

if [ ! -d ${analysis}/deCODE ]; then mkdir -p ${analysis}/deCODE; fi
if [ ! -d ${analysis}/UKB_PPP ]; then mkdir -p ${analysis}/UKB_PPP; fi

module load ceuadmin/R

function caprion_caprion_dr()
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
       z[[i]] <- pQTLtools::novelty_check(k,dr,ldops=list(bfile=bfile,plink=plink),flanking=500000)
    }
    z[["23"]] <- dplyr::mutate(z[["X"]],known.seqnames="23",query.seqnames="23")
    r <- dplyr::filter(dplyr::bind_rows(z[-which(names(z)=="X")]),r2>=0.8) %>%
         dplyr::mutate(seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(panel,caprion,caprion_dr,r,file=file.path(analysis,"reports","caprion-caprion_dr.rda"))
    write.table(r,file=file.path(analysis,"reports","caprion-caprion_dr.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
    cvt <- "~/Caprion/analysis/work/caprion_dr.cis.vs.trans"
    f <- read.csv(cvt) %>%
         mutate(seqname=paste0("chr",SNPChrom),start=as.integer(SNPPos),end=SNPPos) %>%
         arrange(SNPChrom,SNPPos)
END
}

function dr_peptide()
{
  cat <(cat  ${analysis}/peptide/*/fp/*-tbl.tsv | head -1) \
      <(ls ${analysis}/peptide/*/fp/*-tbl.tsv | parallel -C' ' 'sed "1d" {}') > ${analysis}/reports/peptide-tbl.tsv
  R --no-save -q <<\ \ END
    suppressMessages(library(dplyr))
    suppressMessages(library(gap))
    require(stringr)
    options(width=2000)
    analysis <- Sys.getenv("analysis")
    INF <- Sys.getenv("INF")
    plink <- Sys.getenv("plink")
    freq <- read.table("~/Caprion/analysis/bgen/caprion.freq",col.names=c("SNP","REF","ALT","ALT_FREQS")) %>%
            filter(!duplicated(SNP))
    tbl_dr <- read.delim(file.path(analysis,"work","tbl_dr.tsv")) %>%
              dplyr::select(prot,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,Direction,HetISq,logHetP,N) %>%
              dplyr::mutate(Allele1=toupper(Allele1),Allele2=toupper(Allele2)) %>%
              dplyr::rename(Protein=prot,SNP=MarkerName,EA=Allele1,OA=Allele2,EAF=Freq1,Log10P=log.P.)
    peptide_tbl <- read.delim(file.path(analysis,"reports","peptide-tbl.tsv")) %>%
                   dplyr::select(prot,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,Direction,HetISq,logHetP,N) %>%
                   dplyr::mutate(Allele1=toupper(Allele1),Allele2=toupper(Allele2)) %>%
                   dplyr::rename(isotope=prot,SNP=MarkerName,EA=Allele1,OA=Allele2,EAF=Freq1,Log10P=log.P.)
    Protein_dr <- read.csv(file.path(analysis,"work","caprion_dr.cis.vs.trans")) %>%
                  dplyr::mutate(geneRegion=paste(geneChrom,":",geneStart,"-",geneEnd)) %>%
                  dplyr::left_join(freq) %>%
                  dplyr::arrange(prot,SNPChrom,SNPPos) %>%
                  dplyr::rename(Protein=prot) %>%
                  dplyr::select(-log10p,-geneChrom,-geneStart,-geneEnd,-SNPChrom,-SNPPos,-cis) %>%
                  dplyr::left_join(tbl_dr) %>%
                  dplyr::select(Gene,Protein,geneRegion,SNP,Type,EA,OA,EAF,Effect,StdErr,Log10P,Direction,HetISq,logHetP,N,REF,ALT,ALT_FREQS)
    peptide <- read.csv(file.path(analysis,"reports","peptide.cis.vs.trans")) %>%
               dplyr::mutate(geneRegion=paste(geneChrom,":",geneStart,"-",geneEnd)) %>%
               dplyr::arrange(prot,SNPChrom,SNPPos) %>%
               dplyr::rename(Protein=prot) %>%
               dplyr::left_join(peptide_tbl) %>%
               dplyr::left_join(freq) %>%
               dplyr::select(-log10p,-geneChrom,-geneStart,-geneEnd,-SNPChrom,-SNPPos,-cis) %>%
               dplyr::select(Gene,Protein,geneRegion,SNP,Type,EA,OA,EAF,Effect,StdErr,Log10P,Direction,HetISq,logHetP,N,REF,ALT,ALT_FREQS)
    panel <- pQTLdata::caprion %>%
             transmute(prot=gsub("_HUMAN","",Protein),uniprot=Accession)
    dr <- dplyr::mutate(Protein_dr, ProteinSNP=paste0(Protein,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
          dplyr::rename(prot=Protein,rsid=SNP) %>%
          dplyr::left_join(panel)
    peptide <- dplyr::mutate(peptide, ProteinSNP=paste0(peptide,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
               dplyr::rename(prot=Protein,rsid=SNP) %>%
               dplyr::left_join(panel)
    z <- list()
    for(i in unique(dplyr::pull(caprion,chr)))
    {
       j <- dplyr::filter(dr,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       k <- dplyr::filter(peptide,chr %in% i) %>%
            dplyr::select(chr,pos,uniprot,rsid,prot)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       z[[i]] <- pQTLtools::novelty_check(j,k,ldops=list(bfile=bfile,plink=plink),flanking=500000)
    }
    z[["23"]] <- dplyr::mutate(z[["X"]],known.seqnames="23",query.seqnames="23")
    r <- dplyr::filter(dplyr::bind_rows(z[-which(names(z)=="X")]),r2>=0.8) %>%
         dplyr::mutate(seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos) %>%
         unique()
    r %>% mutate(id=paste0(known.prot,"-",known.rsid)) |> select(id) |> unique() |> dim()
    r %>% mutate(id=paste0(query.prot,"-",query.rsid)) |> select(id) |> unique() |> dim()
    r %>%
    group_by(known.prot) %>%
    summarize(n=n(),protein=paste(known.rsid,collapse=","),peptide=paste(query.rsid,collapse=","),r2=paste(r2,collapse=","))
  # .lst files generated from utils.sh
    ZWK_protein1 <- scan(file.path(analysis,"reports","ZWK1.lst"),what="")
    ZYQ_protein1 <- scan(file.path(analysis,"reports","ZYQ1.lst"),what="")
    UDP_protein1 <- scan(file.path(analysis,"reports","UDP1.lst"),what="")
    filter(r,known.prot%in%c(ZWK_protein1,ZYQ_protein1,UDP_protein1)) %>% select(known.prot) |> table() |> length()
    save(panel,dr,peptide,r,file=file.path(analysis,"reports","dr-peptide.rda"))
    write.table(r,file=file.path(analysis,"reports","dr_peptide.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
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
    sentinels <- openxlsx::read.xlsx(xlsx,sheet=2,startRow=3,colNames=TRUE) %>%
               dplyr::select(10,9,3,7,11,14,15,18,19,27,29) %>%
               setNames(c("SeqId","uniprot","gene","prot","variant","chr","pos","A1","A2","beta","log10p")) %>%
               dplyr::filter(A1!="!" & A2!="!") %>%
               dplyr::mutate(chr=gsub("chr","",chr),
                             se=TwoSampleMR::get_se(beta,10^(-log10p))) %>%
               dplyr::mutate(seqname=chr,start=as.integer(pos),end=pos) %>%
               dplyr::arrange(chr,pos)
    gr <- with(sentinels,GenomicRanges::GRanges(seqnames=seqname,
                                                IRanges::IRanges(start,end,names=variant),
                                                rsid=variant,chr38=chr,pos38=pos,A1=A1,A2=A2,
                                                prot=prot,uniprot=uniprot,
                                                beta=beta,se=se,log10p=log10p))
    uniprots <- intersect(with(sentinels,uniprot),caprion_dr[["uniprot"]])
    print(length(uniprots))
    f <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- rtracklayer::import.chain(f)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    deCODE <- rtracklayer::liftOver(gr,ch) %>%
              data.frame %>%
              dplyr::filter(uniprot %in% uniprots) %>%
              dplyr::mutate(seqnames=gsub("chr","",seqnames)) %>%
              dplyr::select(-group_name,-group,-end,-width,-strand) %>%
              dplyr::rename(chr=seqnames,pos=start) %>%
              dplyr::mutate(snpid=gap::chr_pos_a1_a2(chr,pos,A1,A2,prefix=""))
    z <- list()
    for(i in unique(dplyr::pull(caprion_dr,chr)))
    {
       decode <- dplyr::filter(deCODE,chr %in% i) %>%
                 dplyr::select(chr,pos,uniprot,snpid,prot) %>%
                 dplyr::rename(rsid=snpid)
       dr <- dplyr::filter(caprion_dr,chr %in% i) %>%
             dplyr::select(chr,pos,uniprot,rsid,prot)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       z[[i]] <- pQTLtools::novelty_check(decode,dr,ldops=list(bfile=bfile,plink=plink),flanking=500000)
    }
    z[["23"]] <- dplyr::mutate(z[["X"]],known.seqnames="23",query.seqnames="23")
    s <- dplyr::bind_rows(z[-which(names(z)=="X")]) %>%
         dplyr::mutate(seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(deCODE,s,file=file.path(analysis,"deCODE","deCODE.rda"))
    write.table(dplyr::filter(s,r2>=0.8),file=file.path(analysis,"deCODE","deCODE.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
    f <- file.path(Sys.getenv("analysis"),"deCODE","SomaLogicv4.tsv")
    SomaLogicv4 <- openxlsx::read.xlsx(xlsx,sheet=1,startRow=3,colNames=TRUE,cols=1:12)
    write.table(SomaLogicv4,file=f,quote=FALSE,row.names=FALSE,sep="\t")
END
}

function UKB_PPP()
{
  R --no-save -q <<END
    options(width=200)
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
    ST10 <- openxlsx::read.xlsx(url,"ST10",startRow=4) %>%
            dplyr::rename(snpid37=1) %>%
            dplyr::mutate(uniprot=Target.UniProt,rsid=rsID,prot=Assay.Target) %>%
            dplyr::mutate(prot_rsid=paste0(uniprot,"-",rsid))
    uniprots <- intersect(with(ST10,uniprot),with(caprion_dr,uniprot))
    print(length(uniprots))
    UKB_PPP <- dplyr::filter(ST10,uniprot %in% uniprots) %>%
               dplyr::mutate(chrpos=strsplit(snpid37,":"),
                        chr=unlist(lapply(chrpos,"[[",1)),
                        pos=as.integer(unlist(lapply(chrpos,"[[",2))),
                        a1=unlist(lapply(chrpos,"[[",3)),
                        a2=unlist(lapply(chrpos,"[[",4))) %>%
               dplyr::filter(!is.na(a1)&!is.na(a2)) %>%
               dplyr::mutate(rsid=gap::chr_pos_a1_a2(chr,pos,a1,a2,prefix="")) %>%
               dplyr::select(-chrpos)
    z <- list()
    for(i in unique(dplyr::pull(caprion_dr,chr)))
    {
       ukb_ppp <- dplyr::filter(UKB_PPP,chr %in% i) %>%
                  dplyr::select(chr,pos,uniprot,rsid,prot)
       dr <- dplyr::filter(caprion_dr,chr %in% i) %>%
             dplyr::select(chr,pos,uniprot,rsid,prot)
       bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
       z[[i]] <- novelty_check(ukb_ppp,dr,ldops=list(bfile=bfile,plink=plink),flanking=500000)
    }
    z[["23"]] <- dplyr::mutate(z[["X"]],known.seqnames="23",query.seqnames="23")
    t <- dplyr::bind_rows(z[-which(names(z)=="X")]) %>%
         dplyr::filter(r2>=0.8) %>%
         dplyr::mutate(seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
         dplyr::arrange(seqnames,pos) %>%
         dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(UKB_PPP,t,file=file.path(analysis,"UKB_PPP","UKB_PPP.rda"))
    write.table(t,file=file.path(analysis,"UKB_PPP","UKB_PPP.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
END
}

caprion_caprion_dr
deCODE
UKB_PPP
