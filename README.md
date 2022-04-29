 Analysis

This is now done in a named sequence.

* 1_pca_projection.sh
* 2_ggm.R
* 3_wgcna.R
* 4_pca_clustering.R
* 5_pgwas.sh
* 6_meta_analysis.sh
* 7_merge.sh
* 8_hla.R

## 1. Data handling and PCA projection

The pipeline follows HGI contributions nevertheless only serves for reassurance since the study samples were carefully selected. 

## 2. GGM experiment

The results are ready to report.

## 3. WGCNA experiment

This can be finalised according to the Science paper.

## 4. PCA and clustering

The groupings based on proteins can be made on three phases altogether and instead of a classification indicator the first three PCs are used.

The PLINK2 has been consistent in the pilot studies, so `scale()` operaions can be used in the inverse normal transformations.

The phenotypic data is generated in accordance with the double transformations as in SCALLOP-Seq analysis.

## 5. pGWAS

The bgen files were extracted from a list of all samples, the variant IDs of which were replaced when RSid is missing (.).

The bgen generation is moved into .sb based on cclake but can be switched back to cardio by uncommenting the ##SBATCH lines. Also note that

* GCTA/fastGWA employs MAF>=0.0001 (~56%) and geno=0.1 so potentially we can have .bgen files as such to speed up.
* GCTA uses headerless phenotype files, so **the following section from `5_pgwas.sh` is run** in preparation.

    >
    > ```bash
    > sed -i '1d' ${caprion}/work/caprion-1.pheno
    > sed -i '1d' ${caprion}/work/caprion-2.pheno
    > sed -i '1d' ${caprion}/work/caprion-3.pheno
    > ```
    >
    > at `pilot/work` while the original version is saved at `analysis/work/`.

It looked to take 3.5 days on Cardio without unfiltered genotypes but ~12 hours on cclake, and once these are taken care of the analysis can be propagated.

The (sb)atch file is extended to produce Q-Q/Manhattan/LocusZoom plots.

## 6. Meta-analysis

This follows from the SCALLOP/INF implementation, as designed in the logic of a Makefile, i.e.,

```bash
6_meta_analysis <task>
```
where task = METAL_list, METAL_files, METAL_analysis, respectively in sequence.

To extract significant variants one may resort to `awk 'NR==1||$12<log(1e-6)/log(10)' 1433B-1.tbl`, say.

## 7. Variant identification

An iterative merging scheme is employed; the HLA region is simplified but will be specifically handled.

## 8. HLA imputation

This will be based on information from several software, e.g., [CookHLA](https://cambridge-ceu.github.io/csd3/applications/CookHLA.html).
