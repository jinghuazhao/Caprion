# Summary statistics

options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(arrayQualityMetrics))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(VennDiagram))

load("ZWK.rda")
load("ZYQ.rda")
load("UDP.rda")

check <- function(suffix)
{
  cat(suffix,",",sep="")
  protein_suffix <- get(paste0("protein_",suffix))
  mapping_suffix <- get(paste0("mapping_",suffix))
  peptide_suffix <- get(paste0("peptide_",suffix))
  cat("Protein:",length(featureNames(protein_suffix)),",",sep="")
  cat("Protein in mapping:",nrow(subset(mapping_suffix,featureNames(protein_suffix)%in%Protein)),",",sep="")
  cat("Peptide:",length(featureNames(peptide_suffix)),",",sep="")
  cat("Peptide in mapping:",nrow(subset(mapping_suffix,featureNames(peptide_suffix)%in%Isotope.Group.ID)),"\n")
}

check("ZWK")
check("ZYQ")
check("UDP")

# outliers

list_outliers <- function(es, method="upperquartile") outliers(exprs(es),method=method)

for (method in c("KS","sum","upperquartile"))
{
  ZWK_outliers <- list_outliers(protein_ZWK,method=method)
  print(ZWK_outliers@statistic[ZWK_outliers@which])
}

# sumstats
sumstats <- function(es,FUN=median) as.data.frame(t(apply(exprs(es),1,FUN)))
stat_test <- function(FUN)
{
  types <- c("protein","peptide")
  prot_pept <- list()
  for (protpept in 1:2)
  {
    stat_ZWK <- sumstats(get(paste0(types[protpept],"_ZWK")),FUN)
    stat_ZYQ <- sumstats(get(paste0(types[protpept],"_ZYQ")),FUN)
    stat_UDP <- sumstats(get(paste0(types[protpept],"_UDP")),FUN)
    d_a <- data.frame(stat=t(stat_ZWK),batch="pilot")
    d_b <- data.frame(stat=t(stat_ZYQ),batch="batch2")
    d_c <- data.frame(stat=t(stat_UDP),batch="batch3")
    d_long <- bind_rows(d_a,d_b,d_c)
    mm <- model.matrix(~batch,data=d_long)
    m <- lm(stat ~ batch, data=d_long)
    print(summary(m))
    print(car::Anova(m))
    d_long <- mutate(d_long, batch=recode_factor(batch, '0' = 'pilot', '1' = 'batch2', '2' = 'batch3'))
    ggplot(d_long, aes(x = stat)) + theme_bw() + cowplot::theme_cowplot() +
                                    geom_histogram(fill = "white", colour = "black") + facet_grid(batch ~ ., scales = "free") +
                                    ggtitle(paste0(types[protpept],"-",FUN)) + xlab(sub("p","P",types[protpept])) + ylab("Frequency")
    ggsave(paste0(types[protpept],"-",FUN,".pdf"),device="pdf")
    d_wide <- bind_rows(stat_ZWK,stat_ZYQ,stat_UDP)
    rownames(d_wide) <- c("pilot","batch2","batch3")
    pairs(t(d_wide),pch=19)
    title(paste(types[protpept],"(",FUN,")"))
    prot_pept[[protpept]] <- list(long=d_long,wide=d_wide)
  }
  prot_pept
}
pdf("pairs.pdf")
par(mfrow=c(1,2))
means <- stat_test("mean")
sds <- stat_test("sd")
medians <- stat_test("median")
dev.off()

# overlap

## by protein and peptide
overlap <- function(A,B) unlist(lapply(calculate.overlap(list(featureNames(A),featureNames(B))),length))

overlap(protein_ZWK,protein_ZYQ)
overlap(protein_ZWK,protein_UDP)
overlap(protein_ZYQ,protein_UDP)
overlap(peptide_ZWK,peptide_ZYQ)
overlap(peptide_ZWK,peptide_UDP)
overlap(peptide_ZYQ,peptide_UDP)

protein <- list(pilot=featureNames(protein_ZWK),batch2=featureNames(protein_ZYQ),batch3=featureNames(protein_UDP))
dr <- list(pilot=featureNames(dr_ZWK),batch2=featureNames(dr_ZYQ),batch3=featureNames(dr_UDP))
peptide <- list(pilot=featureNames(peptide_ZWK),batch2=featureNames(peptide_ZYQ),batch3=featureNames(peptide_UDP))
unlist(lapply(calculate.overlap(protein),length))
unlist(lapply(calculate.overlap(peptide),length))
VennDiagram::venn.diagram(protein,"protein_ZWK-ZQY-UDP.png",disable.logging=TRUE,height=5,width=5,units="in")
VennDiagram::venn.diagram(peptide,"peptide_ZWK-ZQY-UDP.png",disable.logging=TRUE,height=5,width=5,units="in")
unlist(lapply(calculate.overlap(dr),length))
VennDiagram::venn.diagram(dr,"dr_ZWK-ZQY-UDP.png",disable.logging=TRUE,height=5,width=5,units="in")

## by protein_peptide

protein_peptide <- function(suffix)
{
  mapping <- get(paste0("mapping_",suffix))
  protein_peptide <- select(mapping,Protein,Isotope.Group.ID,Modified.Peptide.Sequence) %>%
                     mutate(protein_peptide=paste0(gsub("_HUMAN","",Protein),"-",Isotope.Group.ID))
}

protein_peptide_ZWK <- protein_peptide("ZWK")
protein_peptide_ZYQ <- protein_peptide("ZYQ")
protein_peptide_UDP <- protein_peptide("UDP")

overlap_protein_peptide <- function(A,B) unlist(lapply(calculate.overlap(list(A$protein_peptide,B$protein_peptide)),length))

overlap_protein_peptide(protein_peptide_ZWK,protein_peptide_ZYQ)
overlap_protein_peptide(protein_peptide_ZWK,protein_peptide_UDP)
overlap_protein_peptide(protein_peptide_ZYQ,protein_peptide_UDP)

protein_peptide_all <- list(pilot=filter(protein_peptide_ZWK,Protein%in%featureNames(protein_ZWK))$protein_peptide,
                            batch2=protein_peptide_ZYQ$protein_peptide,
                            batch3=protein_peptide_UDP$protein_peptide)
unlist(lapply(calculate.overlap(protein_peptide_all),length))
VennDiagram::venn.diagram(protein_peptide_all,"protein_peptide_ZWK-ZYQ-UDP.png",cat.pos=c(0,45,-45),disable.logging=TRUE,height=5,width=5,units="in")

protein_peptide_12 <- right_join(subset(protein_peptide_ZWK,Protein!="-"),protein_peptide_ZYQ,by="protein_peptide")
protein_peptide_23 <- right_join(subset(protein_peptide_ZYQ,Protein!="-"),protein_peptide_UDP,by="protein_peptide")
mapping_12 <- right_join(mapping_ZWK[,1:3],mapping_ZYQ,by="Isotope.Group.ID")
## correct report
dim(subset(protein_peptide_12,is.na(Protein.x)))
## wrong report
dim(subset(protein_peptide_12,Protein.x!=Protein.y))

