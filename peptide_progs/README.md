# Peptide analysis

## CSD3 directory

**/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/**

## Scripts and results

The project directory above contains scripts at **peptide_progs/** and results results at **peptide/**, respectively.

Script name| Description          | Protein-specific error/output
-----------|----------------------|-----------------------------------------------------------
0_utils.sh | Code snippets
1_pgwas.sh | Association analysis[^association] | {protein}.e / {protein}.o
2_meta_analysis.sh | Meta-analysis| {protein}-METAL\_{SLURM\_job\_id}\_{phenotype\_number}.e / {protein}-METAL\_{SLURM\_job\_id}\_{phenotype\_number}.o
3 | Signal identification[^location]
3.1_extract.sh | Signal extraction | \_step1\_{SLURM\_job\_id}\_{phenotype\_number}.e / \_step1\_{SLURM\_job\_id}\_{phenotype\_number}.o
3.2_collect.sh | Signal collection/classification | \_step2\_{protein}.e / \_step2\_{protein}.o
3.3_plot.sh | Forest, Q-Q, Manhattan, LocusZoom, mean-by-genotype/dosage plots | \_step3\_{SLURM\_job\_id}\_{phenotype\_number}.e / \_step3\_{SLURM\_job\_id}\_{phenotype\_number}.o

Prerequistes for a Manhattan/peptide association plot are

- a call to `gz()` (in `0_utils.sh` for protein) for a compressed DR-filtered data.
- Ensembl-VEP (step 3.2 above) but `ceuadmin/ensembl-vep/104` has to be `ceuadmin/ensembl-vep/111-icelake`(so feasible with `vep -i variants.vcf --plugin LoFtool,scores_file.txt`) but plugins have to be added.

The CSD3 icelake module `ceuadmin/R/4.3.3-icelake` now works as smoothly as `ceuadmin/R` for cclake.

## Footnotes

[^association]: **Association**

    Accidentally, they were removed for the benchmarks: A1BG, APOB EPCR, ERAP2, PROC.

[^location]: **Location**

    **{protein}/sentinels/slurm**