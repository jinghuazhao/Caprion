# Peptide analysis

Peptide-level analysis

This mirrors protein-level analysis as in peptide_progs/,

Name       | Description
-----------|----------------------
0_utils.sh | Code snippets
1_pgwas.sh | Association analysis.
2_meta_analysis.sh | Meta-analysis.
3_merge.sh | Signal detection/classification, forest, Q-Q, Manhattan, LocusZoom, mean-by-genotype/dosage plots.
-----------|----------------------

NB. `3_merge.sh` itself is comprised of steps:

#### 3.1. signal extraction
#### 3.2. collection
#### 3.3. graphical representation

which is operationally done through commenting/uncommenting calls to each step exclusively. In particular, CO3 and ITIH2 are resumed after the 12hr threshold for SLURM.

Prerequistes for a Manhattan/peptide association plot by `0_utils.sh` are

- a call to `gz()` for a compressed DR-filtered data; this is to be followed by execution of the R script inside.
- also step 3.2 above requires Ensembl-VEP but `ceuadmin/ensembl-vep/104` does not work on icelake; a remedy is made with `ceuadmin/ensembl-vep/111-icelake` however the loftee plugin is likely to require additional work -- as documented this is feasible with `vep -i variants.vcf --plugin LoFtool,scores_file.txt`.
