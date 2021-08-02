# 27-1-2021 JHZ

uniprot <- Sys.getenv("uniprot")
protein <- Sys.getenv("protein")
gz <- gzfile(paste0("bgen/",uniprot,"-plink2.gz"));
require(qqman);
tbl <- read.delim(gz,as.is=TRUE);
tbl <- within(tbl,{
   SNP <- ID
   CHR <- as.numeric(X.CHROM)
   BP <- POS
})
tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
par(cex=2.5)
qq <- paste0(uniprot,"_qq.png");
png(qq,width=12,height=10,units="in",pointsize=4,res=300)
qq(with(tbl,P),cex.lab=2.5)
dev.off()
manhattan <- paste0(uniprot,"_manhattan.png");
png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
manhattan(tbl,main=uniprot,genomewideline=-log10(8.210181e-12),suggestiveline=FALSE,ylim=c(0,9),cex.axis=2.5,cex.lab=2.5);
dev.off();
