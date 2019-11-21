# 21-11-2019 JHZ

affy <- read.delim("affymetrix.tsv",as.is=TRUE)
ps <- with(affy, phenoscanner::phenoscanner(snpquery=rsid, catalogue="pQTL"))
snps <- with(ps,snps)
results <- with(ps,results)
snps <- within(snps,
{
  chrpos <- sapply(levels(hg19_coordinates),strsplit,":")
  chrom <- unlist(lapply(chrpos,"[",1))
  pos_hg19 <- as.integer(unlist(lapply(chrpos,"[",2)))
  chromStart <- pos_hg19-1
  chromEnd <- pos_hg19
})
results <- within(results,
{
  chrpos <- strsplit(hg19_coordinates,":")
  chrom <- unlist(lapply(chrpos,"[",1))
  pos_hg19 <- as.integer(unlist(lapply(chrpos,"[",2)))
  chromStart <- pos_hg19-1
  chromEnd <- pos_hg19
})
snps_vars <- c("chrom","chromStart","chromEnd","snp","a1","a2","eur","ensembl","hgnc")
result_vars <- c("chrom","chromStart","chromEnd","pmid","beta","se","p","n","efo","snp","trait")
write.table(snps[snps_vars],file="affymetrix.snps",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(results[result_vars],file="affymetrix.results",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
