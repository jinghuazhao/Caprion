# 4-3-2020 JHZ

uniprot <- Sys.getenv("uniprot");
print(uniprot);
gz <- gzfile(paste0("bgen/",uniprot,"-plink2.gz"));
require(qqman);
tbl <- read.delim(gz,as.is=TRUE);
tbl <- within(tbl,{
   SNP <- ID
   CHR <- as.numeric(X.CHROM)
   BP <- POS
})
tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
# qq <- paste0(uniprot,"_qq.png");
# png(qq,width=12,height=10,units="in",pointsize=4,res=300)
# qq(with(tbl,P))
# dev.off()
manhattan <- paste0(uniprot,"_manhattan.png");
png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
manhattan(tbl,main=uniprot,genomewideline=-log10(8.210181e-12),suggestiveline=FALSE,ylim=c(0,25));
dev.off();
