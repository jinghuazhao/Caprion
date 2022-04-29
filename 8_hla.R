
test <- function()
{
  suppressMessages(library(HIBAG))
  hlaLociInfo()
  region <- 500   # kb
  for (hla.id in c("A","B","C","DQA1","DQB1","DRB1"))
  {
    H1 <- HLA_Type_Table[, paste0(hla.id, ".1")]
    H2 <- HLA_Type_Table[, paste0(hla.id, ".2")]
    hla <- hlaAllele(HLA_Type_Table$sample.id, H1, H2, locus=hla.id, assembly="hg19")
    hlatab <- hlaSplitAllele(hla, train.prop=0.3)
    snpid <- hlaFlankingSNP(interval.geno$snp.id, interval.geno$snp.position, hla.id, region*1000, assembly="hg19")
    print(length(snpid))
    id <- hla.id
    assign(id,hla)
  }
# model.fn <- system.file("extdata" ,"ModelList.RData", package="HIBAG")
# OutOfBag.fn <- system.file("extdata" ,"OutOfBag.RData", package="HIBAG")
# model <- get(load(model.fn))
# OutOfBag <- get(load(OutOfBag.fn))
# rm(modellist,mobj)
}
# removing lines containing NA's
# l <- apply(!apply(HLA_Type_Table, 1, is.na),2,all)
# HLA_Type_Table <- HLA_Type_Table[l,]

interval <- function()
{
  bed.fn <- file.path("/home/jhz22/Caprion", "pilot", "data", "merged_imputation.bed")
  fam.fn <- file.path("/home/jhz22/Caprion", "pilot", "data", "merged_imputation.fam")
  bim.fn <- file.path("/home/jhz22/Caprion", "pilot", "data", "merged_imputation.bim")
  hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19",rm.invalid.allele=TRUE, import.chr="6")
}

interval.geno <- interval()
set.seed(100)
train.geno <- hlaGenoSubset(interval.geno,
                            snp.sel=match(snpid, interval.geno$snp.id),
                            samp.sel=match(hlatab$training$value$sample.id,
                            interval.geno$sample.id))
test.geno <- hlaGenoSubset(interval.geno, samp.sel=match(hlatab$validation$value$sample.id, interval.geno$sample.id))
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4, verbose.detail=TRUE)
summary(model)

pred <- hlaPredict(model, interval.geno, type="response+dosage")
summary(pred)
