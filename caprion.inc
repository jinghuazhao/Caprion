# 20-11-2019 JHZ

caprion_xlsx <- function()
{
  library(openxlsx)
  # list
  protein_list <- read.xlsx("Caprion_pilot_protein_list.xlsx")
  # workbook
  wb <- "ZWK_EDR_20191002.xlsx"
  Legend <- read.xlsx(wb, sheet = 1, startRow = 3)
  Samples <- read.xlsx(wb, sheet = 2, startRow = 5)
  Annotations <- read.xlsx(wb, sheet = 3, startRow = 5)
  rawIGs <- read.xlsx(wb, sheet = 4, startRow = 5)
  Normalized_Peptides <- read.xlsx(wb, sheet = 5, startRow = 5)
  Protein_All_Peptides <- read.xlsx(wb, sheet = 6, startRow = 5)
  Protein_DR_Filt_Peptides <- read.xlsx(wb, sheet = 7, startRow = 5)
  save(protein_list,Legend,Samples,Annotations,rawIGs,Normalized_Peptides,Protein_All_Peptides,Protein_DR_Filt_Peptides,file="caprion.rda")
}

caprion_inf <- function()
{
  tmp <- read.delim("inf1.tmp",as.is=TRUE,col.names=c("prot","uniprot"))
  olink_inf1 <- read.delim("olink.inf.panel.annot.tsv",as.is=TRUE)[c("target","target.short","uniprot","panel","hgnc_symbol")]
  inf1 <- merge(tmp,olink_inf1,by="uniprot")
  caprion_chk <- merge(inf1,protein_list[c("Accession","Protein","Gene")],by.x="uniprot",by.y="Accession")
  write.csv(caprion_chk,file="caprion_inf1.chk",row.names=FALSE,quote=FALSE)
  source("olink.inc")
  caprion_inf1 <- merge(protein_list[c("Accession","Protein","Gene")],Inflammation[c("Target","UniProt.No.")],by.x="Accession",by.y="UniProt.No.")
  write.csv(caprion_inf1,file="caprion_inf1.csv",row.names=FALSE,quote=FALSE)
  caprion_cvd2 <- merge(protein_list[c("Accession","Protein","Gene")],CVD_II[c("Target","UniProt.No.")],by.x="Accession",by.y="UniProt.No.")
  write.csv(caprion_cvd2,file="caprion_cvd2.csv",row.names=FALSE,quote=FALSE)
  caprion_cvd3 <- merge(protein_list[c("Accession","Protein","Gene")],CVD_III[c("Target","UniProt.No")],by.x="Accession",by.y="UniProt.No")
  write.csv(caprion_cvd3,file="caprion_cvd3.csv",row.names=FALSE,quote=FALSE)
  caprion_neurology <- merge(protein_list[c("Accession","Protein","Gene")],Neurology[c("Target","UniProt.No")],by.x="Accession",by.y="UniProt.No")
  write.csv(caprion_neurology,file="caprion_neurology.csv",row.names=FALSE,quote=FALSE)
  somalogic <- read.delim("SOMALOGIC_Master_Table_160410_1129info.tsv",as.is=TRUE)
  caprion_somalogic <- merge(protein_list[c("Accession","Protein","Gene")],somalogic[c("Target","UniProt")],by.x="Accession",by.y="UniProt")
  write.csv(unique(caprion_somalogic),file="caprion_somalogic.csv",row.names=FALSE,quote=FALSE)
}

plotfun <- function(col) 
{
  d <- df[,col]
  xlab <- "Individual"
  ylab <- colnames(df)[col]
  plot(d, xlab=xlab, ylab = ylab, main=ylab, type = "p", cex=0.6)
  hist(d, xlab=ylab, main="")
  boxplot(d, horizontal=TRUE, cex=0.6)
}

regfun <- function(col)
{
  p <- df[,col]
  l <- lm(p~sex)
  s <- summary(l)
  c <- with(s,coefficients)
  r <- paste(colnames(df)[col],paste(s$coefficients[,4],collapse="\t"),sep="\t")
  cat(r,"\n",append=TRUE, file="sex.tsv", sep="")
  l <- lm(p~age+sex+bmi)
  s <- summary(l)
  c <- with(s,coefficients)
  r <- paste(colnames(df)[col],paste(s$coefficients[,4],collapse="\t"),sep="\t")
  cat(r,"\n",append=TRUE, file="lm.tsv", sep="")
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
  missing <- read.table("merged_imputation.missing",col.names=c("affymetrix_gwasqc_bl","missing"))
  interval <- merge(d2[,-3],missing,by.x="ID_1",by.y="affymetrix_gwasqc_bl")
  samples <- merge(interval,protein,by.x="ID_1",by.y="affymetrix_gwasqc_bl",all=TRUE)
  eigenvec <- read.delim("merged_imputation.eigenvec")
  covariates <- merge(pheno_protein[c("affymetrix_gwasqc_bl","sex","age","bmi")],eigenvec[,-1],by.x="affymetrix_gwasqc_bl",by.y="IID")
  cat("ID_1 ID_2 missing",names(covariates)[-1], pnames,"\n0 0 0 D C C", rep("C",20), "P P P P P P P P\n",file="SomaLogic.sample")
  write.table(merge(merge(interval,covariates,by.x="ID_1",by.y="affymetrix_gwasqc_bl",all.x=TRUE),
                    protein[,-1],by.x="ID_1",by.y="affymetrix_gwasqc_bl",all.x=TRUE),
              file="SomaLogic.sample",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

ae <- function(X,hidden.layers=c(10,2,10),pdf="ae.pdf")
{
   require(ANN2)
   pdf(pdf)
   AE <- autoencoder(X, hidden.layers, loss.type = 'pseudo-huber',
                     activ.functions = c('tanh','linear','tanh'),
                     batch.size = 8, optim.type = 'adam',
                     n.epochs = 1000, val.prop = 0)
 # Plot loss during training
   plot(AE)
 # Make reconstruction and compression plots
   reconstruction_plot(AE, X)
   compression_plot(AE, X)
 # Reconstruct data and show states with highest anomaly scores
   recX <- reconstruct(AE, X)
   sort(recX$anomaly_scores, decreasing = TRUE)
   dev.off()
   save(AE,recX,file="ae.rda")
}

load("caprion.rda")

# phenotypes
names(Samples) <- c("caprion_id","external_id","comment")
phenotypes <- read.delim("interval_caprion_pilot_samples_phenotype_data.tsv",as.is=TRUE)
phenotypes <- within(phenotypes,{
  age <- agepulse
  sex <- sexpulse
  bmi <- wt_bl/ht_bl/ht_bl
  crp <- crp_bl
  transf <- transf_bl
})
pd <- merge(phenotypes[c("caprion_id","affymetrix_gwasqc_bl","sex","age","bmi","crp","transf")],Samples,by="caprion_id")

td <- Protein_All_Peptides
pap <- merge(protein_list[c("Protein","Accession","Gene")], td, by="Protein")
rownames(td) <- gsub("_HUMAN","",td[,1])
t1 <- t(td[,-1])
td <- data.frame(caprion_id=row.names(t1),data.frame(t1))
pheno_protein <- merge(pd,td,by="caprion_id")
df <- pheno_protein[,-(1:8)]
d <- Protein_All_Peptides[,-1]

excl_id <- read.table("11-3.id", as.is=TRUE, header=TRUE)
eleven <- excl_id[1:11,1]
three <- excl_id[12:14,1]
idx11 <- colnames(d)%in%eleven
idx3 <- colnames(d)%in%three
group <- rep(1,ncol(d))
group[idx11] <- 2
group[idx3] <- 3
col.group=c("black","blue","red")

source("makeRLEboxplot.R")
