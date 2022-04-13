# Analysis

This is now done in a named sequence.

* 1_pca_projection.sh
* 2_ggm.R
* 3_wgcna.R
* 4_pca_clustering.R
* 5_pgwas.sh
* 6_meta_analysis.sh

## 1. Data handling and PCA projection

The pipeline follows HGI contributions nevertheless only serves for reassurance since the study samples were carefully selected. 

## 2. GGM experiment

The results are ready to report.

## 3. WGCNA experiment

This can be finalised according to the Science paper.

## 4. PCA and clustering

The groupings based on proteins can be made on three phases altogether. and instead of a classification indicator the first three PCs are used.

The PLINK2 has been consistent in the pilot studies, so `scale()` operaions can be used in the inverse normal transformations.

The phenotypic data is generated in accordance with the double transformations as in SCALLOP-Seq analysis.

## 5. pGWAS

Note that GCTA/fastGWA employs MAF>=0.0001 (~56%) and geno=0.1 so potentially we can have .bgen files as such to speed up.

It looked to take 4.5 days on Cardio without unfiltered genotypes and once these are taken care of the analysis can be propagated.

## 6. Meta-analysis
