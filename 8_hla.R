options(width=200)
analysis <- "~/Caprion/analysis"
suffix <- Sys.getenv("suffix")
hpc_work <- Sys.getenv("HPC_WORK")
region <- 500 # kb
hlaLoci <- c("A","B","C","DPB1","DQA1","DQB1","DRB1")
seed <- 123456
suppressMessages(library(HIBAG))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))

interval <- function()
{
  setwd(file.path(analysis,"work"))
  interval.geno <- hlaBED2Geno("hla.bed", "hla.fam", "hla.bim", assembly="hg19",rm.invalid.allele=TRUE, import.chr="6")
  assign("interval.geno",interval.geno,envir=.GlobalEnv)
  save(interval.geno,file=file.path(analysis,"work","hla.rda"))
  for (hlaId in hlaLoci)
  {
    snpid <- hlaFlankingSNP(interval.geno$snp.id, interval.geno$snp.position, hlaId, region*1000, assembly="hg19")
    print(length(snpid))
  }
}

HIBAG <- function(hlaId,cohort,reference="HapMap")
{
  cat(hlaId,"\n")
  geno <- get(paste0(cohort,".geno"))
  if (reference=="HapMap")
  {
    model <- get(paste0(hlaId,".model"))
    pred <- hlaPredict(model, geno, type="response+dosage")
  } else {
    model <- hlaModelFromObj(model.list[[hlaId]])
    pred <- predict(model, geno, type="response+prob",cl)
  }
  assign(paste0(cohort,".",hlaId),pred,env=.GlobalEnv)
}

lookup <- function(rsid=NULL)
{
  if (!is.null(rsid)) print(grepl(rsid,lapply(model.list,"[[",3)))
}

hlaAssocTestBatch <- function(hlaLocus,batch)
{
  f <- file.path(analysis,"work",paste0("caprion",suffix,"-",batch,".pheno"))
  d <- read.delim(f,check.names=FALSE)
  id <- pull(d,IID)
  p <- select(d,-FID,-IID)
  dst <- hlaLocus
  dst$value <- subset(dst$value,sample.id %in% id)
  dst$postprob <- dst$postprob[,colnames(dst$postprob) %in% id]
  z <- sapply(names(p),function(col)
       {
         cat(col,"\n",sep="")
         y <- p[[col]]
         hlaAssocTest(dst,formula(y~1),model="additive")
       })
}

hla_main <- function()
{
  cl <- makeCluster(8)
  model.list <- get(load(file.path(analysis,"HLA","HIBAG","AffyAxiomUKB-European-HLA4-hg19.RData")))
  lookup("rs2229092")
  load(file.path(analysis,"work","hla.rda"))
  for (hlaId in hlaLoci) HIBAG(hlaId,"interval","UKBAxiom")
  save(interval.A,interval.B,interval.C,interval.DPB1,interval.DQA1,interval.DQB1,interval.DRB1,
       file=file.path(analysis,"work","interval.hla"))
  f <- file.path(analysis,"work",paste0("caprion",suffix,".pheno"))
  pheno <- read.delim(f,check.names=FALSE)
  load(file=file.path(analysis,"work","interval.hla"))
  h <- r <- list()
  for(hlaLocus in hlaLoci)
  {
    for(i in 1:3)
    {
      z <- hlaAssocTestBatch(get(paste("interval",hlaLocus,sep=".")),i)
      f <- file.path(analysis,"work",paste0(paste("interval",i,hlaLocus,sep="-"),".rda"))
      save(z,file=f)
      r[[i]] <- z
    }
    h[[hlaLocus]] <- r
  }
  save(h,file=file.path(analysis,"work","interval-hlaAsoocTest.rda"))
}

# hla_main()
boosting <- function(id="A",impute=FALSE,seed=100)
# interval.A[...], interval.geno
{
  hla.id <- id
  hla <- get(paste("interval",hla.id,.sep=".")
  set.seed(seed)
  if (! impute)
  {
    hlatab <- hlaSplitAllele(hla, train.prop=0.5)
    names(hlatab)
    summary(hlatab$training)
    summary(hlatab$validation)
    region <- 500
    snpid <- hlaFlankingSNP(interval.geno$snp.id,interval.geno$snp.position, hla.id, region*1000, assembly="hg19")
    length(snpid)
    train.geno <- hlaGenoSubset(interval.geno, snp.sel=match(snpid,interval.geno$snp.id),
    samp.sel=match(hlatab$training$value$sample.id,interval.geno$sample.id))
    test.geno <- hlaGenoSubset(interval.geno, samp.sel=match(hlatab$validation$value$sample.id,interval.geno$sample.id))
    set.seed(100)
    model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=100, verbose.detail=TRUE)
    summary(model)
    pred <- hlaPredict(model, test.geno)
    summary(pred)
    mobj <- hlaModelToObj(model)
    save(mobj, file=file.path(analysis,"work",paste0("HIBAG_model_",hla.id,".RData")))
    save(test.geno, file=file.path(analysis,"work",paste0("testgeno_",hla.id,".RData")))
    save(hlatab, file=file.path(analysis,"work",paste0("HLASplit_",hla.id,".RData")))
    hlaClose(model)
  } else {
    mobj <- get(load(paste0("HIBAG_model_",hla.id,".RData")))
    model <- hlaModelFromObj(mobj)
    test.geno <- get(load(paste0("testgeno",hla.id,".RData")))
    hlatab <- get(load(paste0("HLASplit_",hla.id,".RData"))
    pred <- hlaPredict(model, test.geno, type="response")
    summary(pred)
    (comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model, call.threshold=0))
    (comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model, call.threshold=0.5))
  }
}

load(file.path(analysis,"work","hla.rda"))
load(file.path(analysis,"work","interval.hla"))

boosting(id="A")
boosting(id="A",impute=TRUE)

# Notes

bc58 <- function()
{
  bed.fn <- file.path(hpc_work,"CookHLA","example","1958BC.bed")
  fam.fn <- file.path(hpc_work,"CookHLA","example","1958BC.fam")
  bim.fn <- file.path(hpc_work,"CookHLA","example","1958BC.bim")
  bc58.geno <- hlaBED2Geno(bed.fn,fam.fn,bim.fn,assembly="hg19",rm.invalid.allele=TRUE,import.chr="6")
  assign("bc58.geno",bc58.geno,envir=.GlobalEnv)
}

bc58hatk <- function()
{
  bed.fn <- file.path(hpc_work,"HATK","example","wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.bed")
  fam.fn <- file.path(hpc_work,"HATK","example","wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.fam")
  bim.fn <- file.path(hpc_work,"HATK","example","wtccc_filtered_58C_RA.hatk.300+300.chr6.hg18.bim")
  bc58hatk.geno <- hlaBED2Geno(bed.fn,fam.fn,bim.fn,assembly="hg18",rm.invalid.allele=TRUE,import.chr="6")
  assign("bc58hatk.geno",bc58hatk.geno,envir=.GlobalEnv)
}

HapMap_CEU_model <- function()
{
  set.seed(seed)
  for (hlaId in hlaLoci)
  {
    # when removing lines containing NA's
    # l <- apply(!apply(HLA_Type_Table, 1, is.na),2,all)
    # HLA_Type_Table <- HLA_Type_Table[l,]
    H1 <- HLA_Type_Table[, paste0(hlaId, ".1")]
    H2 <- HLA_Type_Table[, paste0(hlaId, ".2")]
    hla <- hlaAllele(HLA_Type_Table$sample.id, H1, H2, locus=hlaId, assembly="hg19")
    hlatab <- hlaSplitAllele(hla, train.prop=0.3)
    id <- hlaId
    assign(id,hla,envir=.GlobalEnv)
    snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,hlaId, region*1000, assembly="hg19")
    train.geno <- hlaGenoSubset(HapMap_CEU_Geno, snp.sel=match(snpid, HapMap_CEU_Geno$snp.id),
                                samp.sel=match(hlatab$training$value$sample.id, HapMap_CEU_Geno$sample.id))
    test.geno <- hlaGenoSubset(HapMap_CEU_Geno, samp.sel=match(hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))
    hlaModel <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4, verbose.detail=TRUE)
    summary(hlaModel)
    assign(paste0(id,".model"),hlaModel,envir=.GlobalEnv)
    pred <- hlaPredict(A.model, HapMap_CEU_Geno, type="response+dosage")
    summary(pred)
    assign(paste0(id,".pred"),pred,envir=.GlobalEnv)
  }
}

tests <- function()
{
  hlaLociInfo()
  hlaId <- "A"
  plot(model.list[[hlaId]])
  bc58()
  HIBAG(hlaId,"bc58")
# fail with wrong assembly
  HIBAG(hlaId,"bc58hatk")
  interval()
  HapMap_CEU_model()
  HIBAG(hlaId,"interval")

  setwd("~/rds/post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped")
  bed <- "merged_imputation.bed"
  fam <- "merged_imputation.fam"
  bim <- "merged_imputation.bim"
  model.fn <- system.file("extdata" ,"ModelList.RData", package="HIBAG")
  model <- get(load(model.fn))
  OutOfBag.fn <- system.file("extdata" ,"OutOfBag.RData", package="HIBAG")
  OutOfBag <- get(load(OutOfBag.fn))
  rm(modellist,mobj)
}
