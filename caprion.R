# 12-12-2019 JHZ

source("caprion.ini")

tromso_sample <- function()
# replication of TROMSO study
{
  d <- tromso_xlsx()
  pp <- names(table(with(d,sheet3["Protein.Name"])))
  overlap <- colnames(t1)[colnames(t1)%in%pp]
  selected <- with(d, sheet3[["Protein.Name"]] %in% overlap)
  s <- with(d, sheet3)[selected,c("Protein.Name","Ensembl.ID","Sentinel.Variant")]
  write.table(s,file="tromso.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
  rsid <- names(table(with(s,Sentinel.Variant)))
  ps <- phenoscanner::phenoscanner(snpquery=rsid, catalogue="pQTL")
  snps <- with(ps,snps)
  sort(as.numeric(levels(with(snps,chr))))
  prot <- pheno_protein[c("caprion_id","affymetrix_gwasqc_bl",overlap)]
  prot <- within(prot,{for(i in names(prot[,-(1:2)])) assign(paste0(i,"_invn"),gap::invnormal(prot[,i]))})
  prot <- prot[,-ncol(prot)]
  id1_id2_missing_covariates_phenotypes <- merge(id1_id2_missing_covariates,prot[,-1],
                                                 by.x="ID_1",by.y="affymetrix_gwasqc_bl")
  snptest_sample(id1_id2_missing_covariates_phenotypes,"tromso.sample",
                 C=c("age","bmi",paste0("PC",1:20)),
                 D="sex",
                 P=names(prot[,-(1:2)]))
}

peptides_sample <- function()
# analysis on peptides
{
  p <- "ERAP2"
  d <- extract_peptide("ERAP2")
  id <- pheno_protein[,1:2]
  peptides <- merge(id,d,by="caprion_id")
  peptides <- within(peptides,{for(i in names(peptides[,-(1:2)])) assign(paste0(i,"_invn"),gap::invnormal(peptides[,i]))})
  peptides <- peptides[,-ncol(peptides)]
  id1_id2_missing_covariates_phenotypes <- merge(id1_id2_missing_covariates,peptides[,-1],
                                                 by.x="ID_1",by.y="affymetrix_gwasqc_bl")
  snptest_sample(id1_id2_missing_covariates_phenotypes,paste0(p,".sample"),
                 C=c("age","bmi",paste0("PC",1:20)),
                 D="sex",
                 P=names(peptides)[-(1:2)])
}

affymetrix <- function()
{
  affymetrix.id <- with(phenotypes,affymetrix_gwasqc_bl)
  write.table(affymetrix.id[!is.na(affymetrix.id)],file="affymetrix.id",row.names=FALSE,col.names=FALSE,quote=FALSE)
  uniprot <- scan("SomaLogic.uniprot",what="")
  SomaLogic <- subset(pap,Accession%in%uniprot)
  d1 <- t(SomaLogic[,-c(1:3)])
  pnames <- colnames(d1) <- with(SomaLogic,Accession)
  d1 <- data.frame(caprion_id=rownames(d1),round(d1,3))
  protein <- merge(phenotypes[c("caprion_id","affymetrix_gwasqc_bl")],d1,by="caprion_id")
  d2 <- read.table("interval.samples",skip=2,col.names=c("ID_1","ID_2","missing"))
  interval <- merge(d2[,-3],missing,by.x="ID_1",by.y="affymetrix_gwasqc_bl")
  samples <- merge(interval,protein,by.x="ID_1",by.y="affymetrix_gwasqc_bl",all=TRUE)
  cat("ID_1 ID_2 missing",names(covariates)[-1], pnames, paste0(pnames,"_invn"),"\n0 0 0 D C C", rep("C",20), 
      "P P P P P P P P P P P P P P P P\n",file="SomaLogic.sample")
  p <- protein[,-1]
  SomaLogic.sample <- merge(merge(interval,covariates,by.x="ID_1",by.y="affymetrix_gwasqc_bl",all.x=TRUE),
                      p,by.x="ID_1",by.y="affymetrix_gwasqc_bl",all.x=TRUE)
  p1 <- SomaLogic.sample[names(p[,-1])]
  SomaLogic.sample <- within(SomaLogic.sample, {for(i in 1:ncol(p1)) assign(paste0(names(p1)[i],"_invn"), gap::invnormal(p1[,i]))})
  write.table(SomaLogic.sample[,-ncol(SomaLogic.sample)], file="SomaLogic.sample",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

# outliers by AE
prot <- pheno_protein[,-(1:9)]
# ae(prot,hidden.layers=c(100,20,30))
r <- ae_caprion(prot,hidden.layers=c(100,20,30))
idr <- cbind(pheno_protein[c("caprion_id","affymetrix_gwasqc_bl")],mse=rowSums(r)/ncol(prot))
ord <- with(idr, order(mse,decreasing=TRUE))
head(idr[ord,],14)
plot(idr[ord,3],cex=0.4)

# UMAP
library(uwot)
pdf("umap.pdf")
df_umap <- umap(df,n_neighbors = 5, learning_rate = 0.5, pca=50, spread = 3)
plot(df_umap,col="black",cex=0.4)
points(df_umap[(1:197)[idx11],],col="blue",cex=0.4)
points(df_umap[(1:197)[idx3],],col="red",cex=0.4)
title(paste("UMAP with n_neighbors=5, learning_rate=0.5, pca=50, spread=3"))
dev.off()

# PCA
ppc <- with(pheno_protein, prcomp(na.omit(pheno_protein[,-(1:8)]), rank=50, scale=TRUE))
write.table(data.frame(with(pheno_protein,caprion_id),with(ppc,x)),file="pca.txt",sep="\t",quote=FALSE,row.names=FALSE)
pdf("pca.pdf")
screeplot(ppc, npcs=20, type="lines", main="PCA screeplot")
idx11 <- with(pheno_protein,caprion_id)%in%eleven
idx3 <- with(pheno_protein,caprion_id)%in%three
with(ppc, {
  plot(x[,1:2], main="PC1 -- PC2", cex=0.5, pch=21)
  points(x[idx11,1],x[idx11,2], col="blue", cex=0.8, pch=10)
  points(x[idx3,1],x[idx3,2], col="red", cex=0.8, pch=16)
  legend("bottom", legend = c("lower detection", "lipemic"), box.lty = 0, cex = 0.8,
         col = c("red", "blue"), horiz = TRUE, inset = c(0, 1), xpd = TRUE, pch = c(10, 16))
})
biplot(ppc,cex=0.1)
title("biplot")
dev.off()

# RLE plot
pdf("rle.pdf", width=50, height=12)
par(mfrow=c(2,1))
makeRLEplot(d, log2.data=FALSE, groups=group, col.group=col.group, cex=0.3, showTitle=TRUE)
makeRLEplot(d, log2.data=FALSE, cex=0.3, showTitle=TRUE, title="Uncoloured relative log expression (RLE) plot")
dev.off()

# overlaps
protein_list <- within(protein_list,{Protein <- gsub("_HUMAN","",Protein)})
caprion_inf()

# correlation
with(pheno_protein,cor(crp,CRP,use="complete.obs"))
with(pheno_protein,cor(transf,TRFE,use="complete.obs"))

# cluster analysis
km <- kmeans(d,3)
with(km, table(cluster))
library(mclust)
mc <- Mclust(d)
summary(mc)
pdf("mc.pdf")
plot(mc, what = "BIC")
dev.off()

# regression
attach(pheno_protein)
options(width=500)
cat("protein","intercept","sex",sep="\t",file="sex.tsv")
cat("\n",append=TRUE,file="sex.tsv")
cat("protein","intercept","sex","age","bmi",sep="\t",file="lm.tsv")
cat("\n",append=TRUE,file="lm.tsv")
sapply(seq(1, ncol(df)), regfun)
detach(pheno_protein)
pdf("qq.pdf")
sex <- read.delim("sex.tsv",as.is=TRUE)
gap::qqunif(with(sex,sex),cex=0.4)
title("QQ plot for sex")
lm <- read.delim("lm.tsv",as.is=TRUE)
gap::qqunif(with(lm,bmi),cex=0.4)
title("QQ plot for BMI")
dev.off()

# EDA plots
df <- t1
rownames(df) <- gsub("ZWK","",rownames(df))
pdf("scatter-histogram-boxwhisker.pdf")
par(mfrow=c(3,1))
sapply(seq(1, ncol(df)), plotfun)
dev.off()

# Box-Whisker plot
pdf("box.pdf", width=50, height=12)
np <- ncol(df)
xtick <- seq(1, np, by=1)
boxplot(df,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(df), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("raw data")
sdf <- scale(df,scale=FALSE)
boxplot(sdf,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(sdf), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("mean-centred data")
ssdf <- scale(df,scale=TRUE)
boxplot(ssdf,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(ssdf), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("mean-centred and scaled data")
dev.off()
