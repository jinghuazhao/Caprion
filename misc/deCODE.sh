#!/usr/bin/bash

export deCODE=~/rds/results/public/proteomics/deCODE
export analysis=~/Caprion/analysis
export v4=SomaLogicv4.tsv

if [ ! -d ${analysis}/deCODE ]; then mkdir -p ${analysis}/deCODE; fi

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
    ch <- import.chain(path)
    seqlevelsStyle(gr) <- "UCSC"
    gr19 <- liftOver(gr,ch)
    write.table(dplyr::select(data.frame(gr19),-group,-end,-width,-strand) %>%
                dplyr::rename(snpid=group_name,chr=seqnames,pos=start),
                file=file.path(analysis,"deCODE","deCODE.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
END
}
