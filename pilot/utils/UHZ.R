# Phenotypic information on UHZ

options(width=200)
caprion <- Sys.getenv("caprion")
caprion <- ifelse(caprion=="",".",caprion)
load("UHZ.rda")

suppressMessages(library(Biobase))
suppressMessages(library(dplyr))
suppressMessages(library(gap))
suppressMessages(library(parallel))

cluster <- function(data, interactive=FALSE)
{
  pca <- prcomp(na.omit(data), rank=10, scale=TRUE)
  pc1pc2 <- with(pca,x)[,1:2]
  rownames(pc1pc2) <- rownames(data)
  library(mclust)
  mc <- Mclust(pc1pc2,G=2)
  with(mc, {
      png(file.path(caprion,"data4","3d.png"), res=300, width=12, height=10, units="in")
      scatterplot3d::scatterplot3d(with(pca,x[,c(2,1,3)]), 
                                   color=c("blue","red")[classification], 
                                   main="Plot of the PC1, PC2 and PC3", pch=16)
      legend("right", legend=levels(as.factor(classification)), col=c("blue", "red"), pch=16)
      dev.off()
  })
  if (interactive()) with(mc, {rgl::plot3d(with(pca,x[,c(2,1,3)]),col=classification)})
  mc
}

dr_protein <- function(eset,out)
{
  vars <- c("caprion_id","sexPulse","agePulse","ethnicPulse","ht_bl","wt_bl","classification","Affymetrix_gwasQC_bl",paste0("PC",1:20))
  rhs <- paste(setdiff(vars,c("caprion_id","ethnicPulse","ht_bl","wt_bl","Affymetrix_gwasQC_bl")),collapse="+")
  mc_classification <- with(cluster(t(exprs(eset))),classification)
  data_omicsmap <- merge(read.csv("pilotsMap_15SEP2021.csv"),
                   read.csv("INTERVAL_OmicsMap_20210915.csv"),by="identifier",all=TRUE) %>%
                   mutate(caprion_id=if_else(caprion_id=="",CAPRION_BO,caprion_id)) %>%
                   full_join(read.csv("INTERVALdata_15SEP2021.csv"),by="identifier") %>%
                   full_join(data.frame(caprion_id=names(mc_classification),classification=mc_classification),by="caprion_id")
# print(subset(data_omicsmap,is.na(Affymetrix_gwasQC_bl))[c("caprion_id","identifier","Affymetrix_gwasQC_bl")],row.names=FALSE)
  pheno <- merge(data_omicsmap,data.frame(caprion_id=sampleNames(eset)),by="caprion_id",all.y=TRUE) %>%
           left_join(read.delim(file.path(caprion,"data","merged_imputation.eigenvec")),by=c("Affymetrix_gwasQC_bl"="X.FID")) %>%
           select(all_of(vars))
  rownames(pheno) <- with(pheno,caprion_id)
  pData(eset) <- pheno
  validObject(eset)
  r <- sapply(featureNames(eset),function(r)
              invnormal(resid(lm(formula(paste(paste0("invnormal(",sub("(^[0-9])","X\\1",r),")"),"~",rhs)),
                        data=eset,na.action=na.exclude))))
  rownames(r) <- sampleNames(eset)
  d <- mutate(data.frame(r),caprion_id=rownames(r)) %>%
       left_join(pheno[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,names(data.frame(r))) %>%
       filter(!is.na(Affymetrix_gwasQC_bl)) %>%
       mutate(caprion_id=Affymetrix_gwasQC_bl)
  names(d) <- c("FID","IID",gsub("HUMAN","invn",featureNames(eset)))
  write.table(d,file=file.path(caprion,"data4",out),quote=FALSE,row.names=FALSE)
  write.table(d["IID"],file=file.path(caprion,"data4",paste0(gsub("pheno","",out),"ind")),quote=FALSE,col.names=FALSE,row.names=FALSE)
  snptest_sample <- right_join(read.table(file.path(caprion,"data","merged_imputation.missing"),col.names=c("IID","missing")),d) %>%
                    rename(ID_1=FID,ID_2=IID)
  gap::snptest_sample(snptest_sample,file.path(caprion,"data4",paste0(gsub("pheno","",out),"sample")),P=names(snptest_sample[,-(1:3)]))
}

dr_protein(protein_UHZ,"protein.pheno")
dr_protein(dr_UHZ,"dr.pheno")

lm_mclapply <- function(row)
{
 tryCatch({
    d <- eset[row]
    y <- featureNames(d)
    f <- paste(paste0("invnormal(",sub("(^[0-9])","X\\1",y),")"),"~",rhs)
    m <- lm(formula(f),data=d,na.action=na.exclude)
    invnormal(resid(m))
 }, error = function(e) {})
}

lm_parLapply <- function(row)
{
  suppressMessages(library(Biobase))
  suppressMessages(library(gap))
  d <- eset[row]
  y <- featureNames(d)
  f <- paste(paste0("invnormal(",sub("(^[0-9])","X\\1",y),")"),"~",rhs)
  m <- lm(formula(f),data=d,na.action=na.exclude)
  invnormal(resid(m))
}

iinvnormal <- function()
{
  cl <- makeCluster(20)
  clusterExport(cl,"eset",envir=parent.frame())
  clusterExport(cl,"rhs",envir=parent.frame())
  r <- parLapply(cl,1:nrow(get("eset",envir=parent.frame())),lm_parLapply)
  stopCluster(cl)
  return(r)
}

peptide <- function(eset,out)
{
  vars <- c("caprion_id","sexPulse","agePulse","ethnicPulse","ht_bl","wt_bl","Affymetrix_gwasQC_bl",paste0("PC",1:20))
  rhs <- paste(setdiff(vars,c("caprion_id","ethnicPulse","ht_bl","wt_bl","Affymetrix_gwasQC_bl")),collapse="+")
  data_omicsmap <- merge(read.csv("pilotsMap_15SEP2021.csv"),read.csv("INTERVAL_OmicsMap_20210915.csv"),by="identifier",all=TRUE) %>%
                   mutate(caprion_id=if_else(caprion_id=="",CAPRION_BO,caprion_id)) %>%
                   full_join(read.csv("INTERVALdata_15SEP2021.csv"),by="identifier") %>%
                   full_join(data.frame(caprion_id=sampleNames(eset)),by="caprion_id")
  data_omicsmap <- data_omicsmap %>%
                   select(-setdiff(setdiff(names(data_omicsmap),vars),"identifier"))
  print(subset(data_omicsmap,is.na(Affymetrix_gwasQC_bl))[c("caprion_id","identifier","Affymetrix_gwasQC_bl")],row.names=FALSE)
  pheno <- merge(data_omicsmap,data.frame(caprion_id=sampleNames(eset)),by="caprion_id",all.y=TRUE) %>%
           left_join(read.delim(file.path(caprion,"data","merged_imputation.eigenvec")),by=c("Affymetrix_gwasQC_bl"="X.FID")) %>%
           select(all_of(vars))
  rownames(pheno) <- with(pheno,caprion_id)
  pData(eset) <- pheno
  validObject(eset)
# only works on an interactive session
# r <- mclapply(1:nrow(eset),lm_mclapply,mc.cores=20)
  r <- iinvnormal()
  r <- data.frame(r)
  names(r) <- featureNames(eset)
  d <- mutate(data.frame(r),caprion_id=rownames(r)) %>%
       left_join(pheno[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,names(data.frame(r))) %>%
       filter(!is.na(Affymetrix_gwasQC_bl)) %>%
       mutate(caprion_id=Affymetrix_gwasQC_bl)
  names(d) <- c("FID","IID",paste0(featureNames(eset),"_invn"))
  write.table(d,file=file.path(caprion,"data4",out),quote=FALSE,row.names=FALSE)
  write.table(d["IID"],file=file.path(caprion,"data4",paste0(gsub("pheno","",out),"ind")),quote=FALSE,col.names=FALSE,row.names=FALSE)
  snptest_sample <- right_join(read.table(file.path(caprion,"data","merged_imputation.missing"),col.names=c("IID","missing")),d) %>%
                    rename(ID_1=FID,ID_2=IID)
  gap::snptest_sample(snptest_sample,file.path(caprion,"data4",paste0(gsub("pheno","",out),"sample")),P=names(snptest_sample[,-(1:3)]))
}

peptide(peptide_UHZ,"peptide.pheno")
