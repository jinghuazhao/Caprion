library(arrayQualityMetrics)
list_outliers <- function(es, method="upperquartile") outliers(exprs(es),method=method)
sumstats <- function(es,FUN=median) as.data.frame(t(apply(exprs(es),1,FUN)))
load("ZWK.rda")
load("ZYQ.rda")
load("UDP.rda")

# outliers

options(width=200)
for (method in c("KS","sum","upperquartile"))
{
  ZWK_outliers <- list_outliers(protein_ZWK,method=method)
  print(ZWK_outliers@statistic[ZWK_outliers@which])
}

# sumstats
library(dplyr)
library(ggplot2)
library(cowplot)
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
library(VennDiagram)
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
VennDiagram::venn.diagram(protein,"protein_ZWK-ZQY-UDP.png")
VennDiagram::venn.diagram(peptide,"peptide_ZWK-ZQY-UDP.png")
unlist(lapply(calculate.overlap(dr),length))
VennDiagram::venn.diagram(dr,"dr_ZWK-ZQY-UDP.png")
