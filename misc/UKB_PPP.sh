function UKB_PPP
{
  Rscript -e '
     options(width=2000)
     suppressMessages(require(dplyr))
     # SCALLOP-INF list
     METAL <- read.delim(file.path(find.package("pQTLtools"),"tests","INF1.METAL")) %>%
              left_join(gap.datasets::inf1[c("prot","gene")]) %>%
              mutate(prot=gene,prot_rsid=paste0(uniprot,"-",rsid),chr=Chromosome,pos=Position)
     # UKB_PPP list
     require(openxlsx)
     results <- "/rds/project/jmmh2/rds-jmmh2-results/public/proteomics"
     url <- file.path(results,"UKB-PPP","doc","sun22.xlsx")
     ST10 <- read.xlsx(url,"ST10",startRow=3) %>%
             mutate(uniprot=Target.UniProt,rsid=rsID,prot=Assay.Target) %>%
             mutate(prot_rsid=paste0(uniprot,"-",rsid))
     sentinels <- left_join(METAL,ST10,by="prot_rsid") %>%
                  select(prot_rsid,cis.trans,rsID) %>%
                  filter(!is.na(rsID))
     inf1 <- c(with(gap.datasets::inf1,uniprot),with(METAL,uniprot)) %>%
             unique()
     overlap <- filter(ST10,uniprot %in% inf1)
     dim(overlap)
     UKB_PPP <- mutate(overlap,
                chrpos=strsplit(overlap[["Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)"]],":"),
                chr=as.integer(unlist(lapply(chrpos,"[[",1))),
                pos=as.integer(unlist(lapply(chrpos,"[[",2))),
                chrpos=paste(chr,pos,sep=":"))
     suppressMessages(require(GenomicRanges))
     b <- novelty_check(UKB_PPP,METAL)
     replication <- filter(b,r2>=0.8)
     INF <- "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/"
     # write.table(replication,file=file.path(INF,"work","UKB-PPP.txt"),
     #             row.names=FALSE,quote=FALSE,sep="\t")
     replication <- read.delim(file.path(find.package("pQTLtools"),"tests","UKB-PPP.txt")) %>%
                    select(known.seqnames,known.rsid,query.rsid,query.prot)
     variant_list <- unique(c(pull(replication,known.rsid),pull(replication,query.rsid)))
     load(file.path(find.package("pQTLtools"),"tests","novel_data.rda"))
     prot_rsid <- with(novel_data,paste0(prot,"-",rsid))
     prot_rsid_repl <- with(replication,paste0(query.prot,"-",query.rsid))
     left <- setdiff(prot_rsid,prot_rsid_repl)
     # local LD reference panel by chromosome
     # r2 <- LDlinkR::LDmatrix(variant_list,pop="CEU",token=Sys.getenv("LDLINK_TOKEN"))
     plink <- "/rds/user/jhz22/hpc-work/bin/plink"
     b <- list()
     for(i in unique(pull(METAL,Chromosome)))
     {
        u <- filter(UKB_PPP,chr %in% i) %>%
             select(chr,pos,uniprot,rsid,prot)
        m <- filter(METAL,Chromosome %in% i) %>%
             select(chr,pos,uniprot,rsid,prot)
        bfile <- file.path(INF,"INTERVAL","per_chr",paste0("chr",i))
        b[[i]] <- novelty_check(u,m,ldops=list(bfile,plink))
     }
     replication2 <- filter(bind_rows(b), r2>=0.8)
     prot_rsid <- with(novel_data %>%
                  left_join(gap.datasets::inf1[c("prot","gene")]),paste0(gene,"-",rsid))
     prot_rsid_repl <- with(replication2,paste0(query.prot,"-",query.rsid))
     novel <- setdiff(prot_rsid,prot_rsid_repl)
  '
}
