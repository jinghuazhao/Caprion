
options(width=200)
suppressMessages(library(HIBAG))
hlaLociInfo()
caprion <- "/home/jhz22/Caprion/"
hpc_work <- Sys.getenv("HPC_WORK")
region <- 500 # kb
hlaLoci <- c("A","B","C","DQA1","DQB1","DRB1")
seed <- 123456

bc58 <- function()
{
  bed.fn <- file.path(hpc_work,"CookHLA","example","1958BC.bed")
  fam.fn <- file.path(hpc_work,"CookHLA","example","1958BC.fam")
  bim.fn <- file.path(hpc_work,"CookHLA","example","1958BC.bim")
  bc58.geno <- hlaBED2Geno(bed.fn,fam.fn,bim.fn,assembly="hg19",rm.invalid.allele=TRUE,import.chr="6")
  assign("bc58.geno",bc58.geno,envir=.GlobalEnv)
}

interval <- function()
{
  bed.fn <- file.path(caprion, "pilot", "data", "merged_imputation.bed")
  fam.fn <- file.path(caprion, "pilot", "data", "merged_imputation.fam")
  bim.fn <- file.path(caprion, "pilot", "data", "merged_imputation.bim")
  interval.geno <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19",rm.invalid.allele=TRUE, import.chr="6")
  assign("interval.geno",interval.geno,envir=.GlobalEnv)
  save(interval.geno,file=file.path(caprion,"analysis","work","hla.rda"))
  for (hlaId in hlaLoci)
  {
    snpid <- hlaFlankingSNP(interval.geno$snp.id, interval.geno$snp.position, hlaId, region*1000, assembly="hg19")
    print(length(snpid))
  }
}

HapMap_CEU_model <- function()
{
  set.seed(seed)
  for (hlaId in hlaLoci)
  {
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
# model.fn <- system.file("extdata" ,"ModelList.RData", package="HIBAG")
# OutOfBag.fn <- system.file("extdata" ,"OutOfBag.RData", package="HIBAG")
# model <- get(load(model.fn))
# OutOfBag <- get(load(OutOfBag.fn))
# rm(modellist,mobj)
# removing lines containing NA's
# l <- apply(!apply(HLA_Type_Table, 1, is.na),2,all)
# HLA_Type_Table <- HLA_Type_Table[l,]
}

HIBAG <- function(hlaId,cohort)
{
  cat(hlaId,"\n")
  model.id <- get(paste0(hlaId,".model"))
  cohort.geno <- get(paste0(cohort,".geno"))
  cohort.pred <- invisible(hlaPredict(model.id, cohort.geno, type = "response+dosage"))
  assign(paste0(cohort,hlaId),cohort.pred,envir=.GlobalEnv)
}

bc58()
interval()
load(file.path(caprion,"analysis","work","hla.rda"))
HapMap_CEU_model()
HIBAG("B","bc58")
HIBAG("B","interval")
