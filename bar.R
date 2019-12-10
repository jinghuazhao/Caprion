# 10-12-2019 JHZ

bar <- function(data, title)
{
  affy <- read.delim(data,as.is=TRUE)
  affy_vars <- c("chromosome","position","add_beta_1","add_se_1","add_pvalue","A","B","uniprot")
  box <- read.table("SomaLogic.box",as.is=TRUE,header=TRUE,fill=TRUE)
  box <- within(box,{Allele1=toupper(Allele1);Allele2=toupper(Allele2)})
  box_vars <- c("chr","pos","Effect","StdErr","logP","protein","SOMAS_ID_round1","rsid","Allele1","Allele2")
  affybox <- merge(box[box_vars],affy[affy_vars],by.x=c("chr","pos"),by.y=c("chromosome","position"))
  affybox <- within(affybox[-5,],{switch <-A!=Allele1;add_beta_1[switch] <- -add_beta_1[switch]})
  SomaLogic <- data.frame(affybox[c("Effect","StdErr","logP","protein","uniprot","rsid")],platform="SomaLogic")
  caprion <- data.frame(affybox[c("add_beta_1","add_se_1","add_pvalue","protein","uniprot","rsid")],platform="Caprion")
  names(caprion)[1:3] <- c("Effect","StdErr","logP")
  caprion <- within(caprion, {logP <- log10(logP)})
  tabbedData <- within(rbind(SomaLogic,caprion),{protein <- paste(protein,rsid,sep="-")})
  limits <- with(tabbedData, aes(ymax = Effect + StdErr, ymin = Effect - StdErr))
  p <- ggplot(data = tabbedData, aes(x = factor(protein), y = Effect, fill = factor(platform)))
  p + geom_bar(stat = "identity", position = position_dodge(0.9)) + geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) + labs(x = "Protein-SNP pair", y = "Effect size") + ggtitle(title) + scale_fill_discrete(name = "platform") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

library(ggplot2)
pdf("bar.pdf")
bar("affymetrix.tsv","Effect size by protein and platform (Caprion untransformed)")
bar("affymetrix_invn.tsv","Effect size by protein and platform (Caprion invnorm)")
dev.off()
