# Peptide analysis

## Scripts / Logs

Peptide-level analysis mirrors protein-level analysis as in **peptide_progs/**,

Name       | Description          | Error/Output per Protein
-----------|----------------------|-----------------------------------------------------------
0_utils.sh | Code snippets
1_pgwas.sh | Association analysis[^association] | {protein}.e/{protein}.o
2_meta_analysis.sh | Meta-analysis| {protein}-METAL\_{SLURM\_job\_id}\_{phenotype\_number}.o
3 | Signal identification[^location]
3.1_extract.sh | Signal extraction | \_step1\_{SLURM\_job\_id}\_{phenotype\_number}.e/\_step1\_{SLURM\_job\_id}\_{phenotype\_number}.o
3.2_collect.sh | Signal collection/classification | \_step2\_{protein}.e/\_step2\_{protein}.o
3.3_plot.sh | Various plots[^plots] | \_step3\_{SLURM\_job\_id}\_{phenotype\_number}.e/\_step1\_{SLURM\_job\_id}\_{phenotype\_number}.o

Prerequistes for a Manhattan/peptide association plot by `0_utils.sh` are

- a call to `gz()` for a compressed DR-filtered data.
- Ensembl-VEP (also step 3.2 above) but `ceuadmin/ensembl-vep/104` has to be `ceuadmin/ensembl-vep/111-icelake` so is feasible with `vep -i variants.vcf --plugin LoFtool,scores_file.txt`.

CO3 and ITIH2 initially required re-run after the 12hr threshold for SLURM, but now both appeared furnished within 12hr.

The CSD3 icelake module `ceuadmin/R/4.3.3-icelake` now works as smoothly as `ceuadmin/R` for cclake

## Footnotes

[^association]: **Association**

    Accidentally, they were removed for the benchmarks: A1BG, APOB EPCR, ERAP2, PROC.

[^location]: **Location**

    **sentinels/slurm**

[^plots]: **Plots**

    These include forest, Q-Q, Manhattan, LocusZoom, mean-by-genotype/dosage plots.
