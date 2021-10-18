# Phenotypic information on UDP

options(width=200)
caprion <- Sys.getenv("caprion")
caprion <- ifelse(caprion=="",".",caprion)
load("UDP.rda")

suppressMessages(library(Biobase))
suppressMessages(library(dplyr))
suppressMessages(library(gap))

cluster <- function(data, interactive=FALSE)
{
  pca <- prcomp(na.omit(data), rank=10, scale=TRUE)
  pc1pc2 <- with(pca,x)[,1:2]
  rownames(pc1pc2) <- rownames(data)
  library(mclust)
  mc <- Mclust(pc1pc2,G=2)
  if (interactive) with(mc,
  {
       png(file.path(caprion,"data3","3d.png"), res=300, width=12, height=10, units="in")
       scatterplot3d::scatterplot3d(with(pca,x[,c(2,1,3)]), color=c("blue","red")[classification], main="Plot of the PC1, PC2 and PC3", pch=16)
       legend("right", legend=levels(as.factor(classification)), col=c("blue", "red"), pch=16)
       dev.off()
       rgl::plot3d(with(pca,x[,c(2,1,3)]),col=classification)
  })
  mc
}

dr_protein <- function(eset,out)
{
  vars <- c("caprion_id","sexPulse","agePulse","ethnicPulse","ht_bl","wt_bl","Affymetrix_gwasQC_bl","classification",paste0("PC",1:20))
  mc_classification <- with(cluster(t(exprs(eset))),classification)
  data_omicsmap <- merge(read.csv("pilotsMap_15SEP2021.csv"),read.csv("INTERVAL_OmicsMap_20210915.csv"),by="identifier",all=TRUE) %>%
                   mutate(caprion_id=if_else(caprion_id=="",CAPRION_BO,caprion_id)) %>%
                   full_join(read.csv("INTERVALdata_15SEP2021.csv"),by="identifier") %>%
                   full_join(data.frame(caprion_id=names(mc_classification),classification=mc_classification),by="caprion_id")
  print(subset(data_omicsmap,is.na(Affymetrix_gwasQC_bl))[c("caprion_id","identifier","Affymetrix_gwasQC_bl")],row.names=FALSE)
  pheno <- merge(data_omicsmap,data.frame(caprion_id=sampleNames(eset)),by="caprion_id",all.y=TRUE) %>%
           left_join(read.delim(file.path(caprion,"data","merged_imputation.eigenvec")),by=c("Affymetrix_gwasQC_bl"="X.FID")) %>%
           select(all_of(vars))
  rownames(pheno) <- with(pheno,caprion_id)
  pData(eset) <- pheno
  validObject(eset)
  rhs <- paste(setdiff(vars,c("caprion_id","ethnicPulse","ht_bl","wt_bl","Affymetrix_gwasQC_bl")),collapse="+")
  r <- sapply(featureNames(eset),function(r) 
              invnormal(resid(lm(formula(paste(paste0("invnormal(",sub("(^[0-9])","X\\1",r),")"),"~",rhs)),data=eset,na.action=na.exclude))))
  rownames(r) <- sampleNames(eset)
  d <- mutate(data.frame(r),caprion_id=rownames(r)) %>%
       left_join(pheno[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,names(data.frame(r))) %>%
       filter(!is.na(Affymetrix_gwasQC_bl)) %>%
       mutate(caprion_id=Affymetrix_gwasQC_bl)
  names(d) <- c("FID","IID",gsub("HUMAN","invn",featureNames(eset)))
  write.table(d,file=file.path(caprion,"data3",out),quote=FALSE,row.names=FALSE,sep="\t")
}

peptide <- function(eset,out)
{
  vars <- c("caprion_id","sexPulse","agePulse","ethnicPulse","ht_bl","wt_bl","Affymetrix_gwasQC_bl",paste0("PC",1:20))
  data_omicsmap <- merge(read.csv("pilotsMap_15SEP2021.csv"),read.csv("INTERVAL_OmicsMap_20210915.csv"),by="identifier",all=TRUE) %>%
                   mutate(caprion_id=if_else(caprion_id=="",CAPRION_BO,caprion_id)) %>%
                   full_join(read.csv("INTERVALdata_15SEP2021.csv"),by="identifier")
  print(subset(data_omicsmap,is.na(Affymetrix_gwasQC_bl))[c("caprion_id","identifier","Affymetrix_gwasQC_bl")],row.names=FALSE)
  pheno <- merge(data_omicsmap,data.frame(caprion_id=sampleNames(eset)),by="caprion_id",all.y=TRUE) %>%
           left_join(read.delim(file.path(caprion,"data","merged_imputation.eigenvec")),by=c("Affymetrix_gwasQC_bl"="X.FID")) %>%
           select(all_of(vars))
  rownames(pheno) <- with(pheno,caprion_id)
  pData(eset) <- pheno
  validObject(eset)
  rhs <- paste(setdiff(vars,c("caprion_id","ethnicPulse","ht_bl","wt_bl","Affymetrix_gwasQC_bl")),collapse="+")
  r <- sapply(featureNames(eset),function(r)
              invnormal(resid(lm(formula(paste(paste0("invnormal(",sub("(^[0-9])","X\\1",r),")"),"~",rhs)),data=eset,na.action=na.exclude))))
  rownames(r) <- sampleNames(eset)
  d <- mutate(data.frame(r),caprion_id=rownames(r)) %>%
       left_join(pheno[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,names(data.frame(r))) %>%
       filter(!is.na(Affymetrix_gwasQC_bl)) %>%
       mutate(caprion_id=Affymetrix_gwasQC_bl)
  names(d) <- c("FID","IID",paste0(featureNames(eset),"_invn"))
  write.table(d,file=file.path(caprion,"data3",out),quote=FALSE,row.names=FALSE,sep="\t")
}

dr_protein(protein_UDP,"protein.tsv")
dr_protein(dr_UDP,"dr.tsv")
peptide(peptide_UDP,"peptide.tsv")

# --- legacy code ---

checks <- function()
{
  head(exprs(protein_UDP))
  head(featureNames(protein_UDP))
  head(sampleNames(protein_UDP))
  experimentData(protein_UDP)
  intersect(OmicsMap$caprion_id,sampleNames(protein_UDP))
  intersect(OmicsMap$Affymetrix_gwasQC_bl,pData(protein_UDP)$Affymetrix_gwasQC_bl)
  nrows <- length(featureNames(protein_UDP))
  ncols <- length(intersect(OmicsMap$caprion_id,sampleNames(protein_UDP)))
}

array_transpose <- function(x)
{
  d <- x[,-1]
  rownames(d) <- x[,1]
  td <- t(d)
}

na_list <- function()
# sed '1d' ${caprion}/data3/protein.tsv | grep -n NA | cut -f1 | sed 's/:NA//' | tr '\n' ' ' > ${caprion}/data3/na_UDP.lst
{
  if(file.exists(file.path(caprion,"data3","na_UDP.lst"))) na_UDP <- scan(file.path(caprion,"data3","na_UDP.lst"))
  exprs(protein_UDP)[,na_UDP]
  t(d[na_UDP,])
}
