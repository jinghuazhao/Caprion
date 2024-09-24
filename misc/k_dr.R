suppressMessages(library(gap))
require(stringr)
options(width=2000)
options("openxlsx.borderColour"="#4F80BD")
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
for(i in unique(pull(caprion,chr)))
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
