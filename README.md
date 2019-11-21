# Caprion pilot study

## Overlaps

* [caprion_inf1.csv](caprion_inf1.csv). Overlap with Olink/Inflammation.
* [caprion_cvd2.csv](caprion_cvd2.csv). Overlap with Olink/CVD II.
* [caprion_cvd3.csv](caprion_cvd3.csv). Overlap with Olink/CVD III.
* [caprion_neurology.csv](caprion_neurology.csv). Overlap with Olink/Neurology.
* [caprion_somalogic.csv](caprion_somalogic.csv). Overlap with SomaLogic.

These are according to the following two files,

* [Olink validation data all panels.xlsx](https://github.com/jinghuazhao/INF/blob/master/doc/Olink%20validation%20data%20all%20panels.xlsx?raw=true).
* [SOMALOGIC_Master_Table_160410_1129info.tsv](https://github.com/jinghuazhao/SomaLogic/blob/master/doc/SOMALOGIC_Master_Table_160410_1129info.tsv?raw=true).

## Correlations

Given CRP and Transferrin from alternative measurements,
```log
> # correlation
> with(pheno_protein,cor(crp,CRP,use="complete.obs"))
[1] 0.3023142
> with(pheno_protein,cor(transf,TRFE,use="complete.obs"))
[1] 0.3616421
```

## Preliminary results

* [box.pdf](box.pdf). Box-Whisker plots for all proteins.
* [rle.pdf](rle.pdf). RLE plot<sup>**1**</sup>.
* [pca.pdf](pca.pdf). PCA analysis
* [scatter-histogram-boxwhisker.pdf](scatter-histogram-boxwhisker.pdf). EDA plots.
* [sex.tsv](sex.tsv). Results (p values) of linear model protein ~ sex.
* [lm.tsv](lm.tsv). Results (p values) of linear model protein ~ sex + age + bmi.
* [qq.pdf](qq.pdf). QQ plots.
* [umap.pdf](umap.pdf). Uniform Manifold Approximation and Projection (UMAP) plot<sup>**2**</sup>.

## Affymetrix data

Files interval.samples and bgen/bgen.bgi are from /DO-NOT-MODIFY-SCRATCH/curated_genetic_data/interval/imputed/.

[affymetrix.tsv](affymetrix.tsv) contains association results ([original comments](SomaLogic.txt)).

The analysis were done using SNPTEST v2.5.5 score method<sup>**3**</sup> for additive model with adjustment for sex, age, bmi, and PC1-PC20 derived from directly genotyped variants from above.

___
**NOTES**

<sup>**1**</sup>Citation: Gandolfo LC, Speed TP (2018) RLE plots: Visualizing unwanted variation in high dimensional data. *PLoS ONE* 13(2): e0191629. https://doi.org/10.1371/journal.pone.0191629

Given protein by individual data (987 x 200 matrix) in object `d`, an RLE plot is done effectively as follows,
```r
medians  <- apply(d, 1, median, na.rm=TRUE)
rle <- sweep(d, MARGIN=1, STATS=medians, FUN='-')
boxplot(rle,xaxt="n",cex=0.2)
xtick <- seq(1, ncol(rle), by=1)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(d), srt=90, pos=1, xpd = TRUE, cex=0.4)
mtext("Relative log expression (RLE) plot",side=3,line=0.1)
```
which is now wrappup as an R function [makeRLEplot.R](makeRLEboxplot.R).

[11-3.id](11-3.id) contains individuals flagged (as *lipemic* [1-11] and *lower detection* [12-14]) by Caprion.

<sup>**2**</sup>Diaz-Papkovich A, Anderson-Trocme L, Ben-Eghan C, Gravel S. UMAP reveals cryptic population structure and phenotype heterogeneity in 
large genomic cohorts. *PLoS Genet*. 2019 Nov 1;15(11):e1008432. doi: 10.1371/journal.pgen.1008432. eCollection 2019 Nov.

<sup>**3**</sup>BOLT-LMM v2.3.4 was able to produced results for the following pairs.

uniprot-protein-SNP | CHR | BP | GENPOS | ALLELE1 | ALLELE0 | A1FREQ | INFO | BETA | SE | P_BOLT_LMM_INF
--------------------|-----|----|--------|---------|---------|--------|------|------|----|---------------
Q6P179-ERAP2-rs2927608 | 5 | 96252432 | 0 | G | A | 0.55102 | 1 | -1.07763 | 0.0971119 | 1.3E-28
P0DJI8-SAA1-rs35179000 | 11 | 18290903 | 0 | T | C | 0.245153 | 0.944861 | 0.0218493 | 0.0386491 | 5.7E-01
