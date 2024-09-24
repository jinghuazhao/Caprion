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
    write.table(select(df,snpid,prot,Pval),file=file.path(INF,"deCODE","replication.tsv"),
                row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    df_ids <- read.table(file.path(INF,"deCODE","olink_inf.full"),
                         col.names=c("uniprot","prot","gene","chrpos","rsid","snpid")) %>%
              left_join(df) %>%
              select(-c(group,rsids,chrpos)) %>%
              select(gene,rsid,chr,pos,effectAllele,otherAllele,ImpMAF,Beta,SE,mlog10p,Pval) %>%
              left_join(select(pQTLdata::inf1,gene,target.short)) %>%
              arrange(gene,rsid,Pval) %>%
              group_by(gene,rsid) %>%
              slice_head(n=1)
    deCODE <- select(df_ids,gene,rsid,chr,pos,effectAllele,otherAllele,Beta,SE,Pval,ImpMAF,mlog10p)
    INF1_METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
                  left_join(select(pQTLdata::inf1,prot,gene)) %>%
                  mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))
    tbl <- select(INF1_METAL,prot,rsid,Chromosome,Position,cis.trans) %>%
           left_join(select(pQTLdata::inf1,gene,prot,target.short)) %>%
           left_join(select(deCODE,gene,rsid,effectAllele,otherAllele,ImpMAF,Beta,SE,Pval)) %>%
           mutate(Protein=target.short,Pval=gap::pvalue(Beta/SE),Pval=if_else(Pval=="NAeNA","NA",Pval)) %>%
           select(Protein,Chromosome,Position,rsid,effectAllele,otherAllele,ImpMAF,Beta,SE,Pval)
    write.table(tbl,file=file.path(INF,"deCODE","deCODE.csv"),row.names=FALSE,quote=FALSE,sep=",")
    INF1_aristotle <- read.delim(file.path(INF,"aristotle","INF1.merge.replication.txt-rsid")) %>%
                      mutate(b_ARISTOTLE=BETA,se_ARISTOTLE=SE,P_ARISTOTLE=PVAL) %>%
                      mutate(prot=unlist(lapply(strsplit(Protein,"-"),"[",1)),
                             rsid=unlist(lapply(strsplit(Protein,"-"),"[",2)),) %>%
                      left_join(select(pQTLdata::inf1,prot,gene,target.short)) %>%
                      left_join(select(INF1_METAL,prot,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans)) %>%
                      select(target.short,gene,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans,
                             EFFECT_ALLELE,REFERENCE_ALLELE,b_ARISTOTLE,se_ARISTOTLE,P_ARISTOTLE)
    INF1_deCODE <- left_join(INF1_METAL,deCODE) %>%
                   mutate(b_deCODE=Beta,se_deCODE=SE,P_deCODE=Pval,mlog10p=-gap::log10p(b_deCODE/se_deCODE)) %>%
                   select(gene,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,effectAllele,otherAllele,-Beta,-SE,-Pval,
                          mlog10p,cis.trans,b_deCODE,se_deCODE,P_deCODE) %>%
                   mutate(sw=if_else(Allele1==effectAllele,1,-1),b_deCODE=sw*b_deCODE,sw2=(sign(Effect)==sign(b_deCODE))+0,
                   col=case_when(
                            mlog10p >= -log10(5e-8) ~ "red",
                            mlog10p >= -log10(0.05/180) & mlog10p <= -log10(5e-8) ~ "orange",
                            mlog10p  < -log10(0.05/180) ~ "grey",
                            TRUE ~ "white" )) %>%
                   filter(!is.na(b_deCODE))
    subset(INF1_deCODE[c("Effect","b_deCODE","log.P.","P_deCODE","cis.trans")],cis.trans=="cis") %>% arrange(Effect)
    nrow(INF1_deCODE)
    filter(INF1_deCODE[c("gene","Effect","b_deCODE","Allele1","Allele2","effectAllele","otherAllele","sw","log.P.","mlog10p")],
           mlog10p>=-log10(5e-8)) %>%
    filter(sign(Effect)!=sign(b_deCODE))
    filter(INF1_deCODE[c("gene","Effect","b_deCODE","Allele1","Allele2","effectAllele","otherAllele","sw","log.P.","mlog10p")],
           mlog10p>=-log10(5e-2/180)) %>%
    filter(sign(Effect)!=sign(b_deCODE))
    INF1_aristotle_deCODE <- left_join(select(INF1_aristotle,target.short,gene,rsid,Allele1,Allele2,Freq1,Effect,StdErr,cis.trans,
                                              EFFECT_ALLELE,REFERENCE_ALLELE,b_ARISTOTLE,se_ARISTOTLE,P_ARISTOTLE),
                                       select(INF1_deCODE,gene,rsid,effectAllele,otherAllele,b_deCODE,se_deCODE,P_deCODE,mlog10p))
    write.table(INF1_aristotle_deCODE,file=file.path(INF,"deCODE","INF1_aristotle_deCODE.csv"),row.names=FALSE,quote=FALSE,sep=",")
    filter(INF1_aristotle, is.na(P_ARISTOTLE)) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 5e-10) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 5e-8) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 0.05/180) %>% nrow()
    filter(INF1_aristotle, P_ARISTOTLE <= 5e-2) %>% nrow()
    filter(INF1_deCODE, is.na(mlog10p)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(5e-10)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(5e-8)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(0.05/180)) %>% nrow()
    filter(INF1_deCODE, mlog10p>=-log10(5e-2)) %>% nrow()
    filter(INF1_aristotle_deCODE,is.na(P_ARISTOTLE) & is.na(mlog10p)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-10 | mlog10p >= -log10(5e-10)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-8 | mlog10p >= -log10(5e-8)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=0.05/180 | mlog10p >= -log10(0.05/180)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-2 | mlog10p >= -log10(5e-2)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE > 5e-8 & mlog10p >= -log10(5e-8)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE <=5e-8 & mlog10p <  -log10(5e-8)) %>% nrow()
    filter(INF1_aristotle_deCODE,P_ARISTOTLE > 5e-8 & mlog10p <  -log10(5e-8)) %>% nrow()
    filter(select(INF1_aristotle_deCODE,target.short,gene,rsid,Allele1,Allele2,effectAllele,otherAllele,
                  b_deCODE,se_deCODE,P_deCODE,P_ARISTOTLE,cis.trans),
           P_ARISTOTLE <=5e-8 & P_deCODE > 5e-8)
    SF_INF_deCODE <- function(d,f)
    {
      png(file=file.path(INF,"deCODE",f),res=300,width=15,height=15,units="in")
      attach(d)
      par(mar=c(5,5,1,1))
      plot(Effect,b_deCODE,pch=19,cex=2,col=col,xaxt="n",ann=FALSE,cex.axis=2)
      axis(1,cex.axis=2,lwd.tick=0.5)
      legend(x=1,y=-1,c("P<=5e-8","P >5e-8 & P <= 0.05/180","P>0.05/180"),
             box.lwd=0,cex=2,col=c("red","orange","grey"),pch=19)
      mtext("deCODE",side=2,line=3,cex=2)
      mtext("Meta-analysis",side=1,line=3,cex=2)
      abline(h=0,v=0)
      detach(d)
      dev.off()
    }
    all <- INF1_deCODE %>%
           filter(!is.na(b_deCODE)) %>%
           select(Effect,b_deCODE)
    cor(all,use="everything")
    cis <- INF1_deCODE %>%
           filter(!is.na(b_deCODE) & cis.trans=="cis") %>%
           select(Effect,b_deCODE)
    cor(cis,use="everything")
    trans <- INF1_deCODE %>%
             filter(!is.na(b_deCODE) & cis.trans=="trans") %>%
             select(Effect,b_deCODE)
    cor(trans,use="everything")
    SF_INF_deCODE(INF1_deCODE,"SF-INF-deCODE.png")
    SF_INF_deCODE(filter(INF1_deCODE,cis.trans=="cis"),"SF-INF-deCODE-cis.png")
    SF_INF_deCODE(filter(INF1_deCODE,cis.trans=="trans"),"SF-INF-deCODE-trans.png")
END
}
