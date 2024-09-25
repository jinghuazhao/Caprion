#!/usr/bin/bash

export deCODE=~/rds/results/public/proteomics/deCODE
export analysis=~/Caprion/analysis
export v4=SomaLogicv4.tsv

if [ ! -d ${analysis}/deCODE ]; then mkdir -p ${analysis}/deCODE; fi
if [ ! -d ${analysis}/UKB_PPP ]; then mkdir -p ${analysis}/UKB_PPP; fi

function caprion()
{
  R --no-save -q <<END
    suppressMessages(library(gap))
    require(stringr)
    options(width=2000)
    analysis <- Sys.getenv("analysis")
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
    INF <- "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF"
    plink <- "/rds/user/jhz22/hpc-work/bin/plink"
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
    b[["23"]] <- mutate(b[["X"]],known.seqnames="23",query.seqnames="23")
    replication <- dplyr::filter(bind_rows(b[-which(names(b)=="X")]),r2>=0.8) %>%
                   dplyr::rename(known.gene=known.uniprot,query.gene=query.uniprot)
    Replication <- dplyr::mutate(replication,seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
                   dplyr::arrange(seqnames,pos) %>%
                   dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
    save(file.path(analysis,"reports","caprion-caprion_dr.rda"))
END
}

function deCODE()
{
  R --no-save -q <<END
    options(width=200)
    suppressMessages(library(dplyr))
    cvt <- "~/Caprion/analysis/work/caprion_dr.cis.vs.trans"
    f <- read.csv(cvt) %>%
         mutate(seqname=paste0("chr",SNPChrom),start=as.integer(SNPPos),end=SNPPos) %>%
         arrange(SNPChrom,SNPPos)
    f <- file.path(Sys.getenv("deCODE"),"doc","ferkingstad21.xlsx")
    SomaLogicv4 <- openxlsx::read.xlsx(f,sheet=1,startRow=3,colNames=TRUE,cols=1:12)
    out <- file.path(Sys.getenv("analysis"),"deCODE",Sys.getenv("v4"))
    write.table(SomaLogicv4,file=out,quote=FALSE,row.names=FALSE,sep="\t")
    sentinels <- openxlsx::read.xlsx(f,sheet=2,startRow=3,colNames=TRUE) %>%
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
    gr19 <- rtracklayer::liftOver(gr,ch)
    write.table(dplyr::select(data.frame(gr19),-group,-end,-width,-strand) %>%
                dplyr::rename(snpid=group_name,chr=seqnames,pos=start),
                file=file.path(analysis,"deCODE","deCODE.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
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
     INF <- "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/"
#    replication <- dplyr::filter(b,r2>=0.8)
     # write.table(replication,file=file.path(INF,"work","UKB-PPP.txt"),
     #             row.names=FALSE,quote=FALSE,sep="\t")
     replication <- read.delim(file.path(find.package("pQTLtools"),"tests","UKB-PPP.txt")) %>%
                    dplyr::select(known.seqnames,known.rsid,query.rsid,query.prot)
     variant_list <- unique(c(dplyr::pull(replication,known.rsid),
                              dplyr::pull(replication,query.rsid)))
     plink <- "/rds/user/jhz22/hpc-work/bin/plink"
     b <- list()
     for(i in unique(dplyr::pull(caprion_dr,Chromosome)))
     {
        u <- dplyr::filter(UKB_PPP,chr %in% i) %>%
             dplyr::select(chr,pos,uniprot,rsid,prot)
        m <- dplyr::filter(METAL,Chromosome %in% i) %>%
             dplyr::select(chr,pos,uniprot,rsid,prot)
        bfile <- file.path(INF,"INTERVAL","per_chr",paste0("chr",i))
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
# deCODE
# UKB_PPP
