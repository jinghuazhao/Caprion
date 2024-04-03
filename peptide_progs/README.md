# Peptide analysis

Peptide-level analysis mirrors protein-level analysis as in **peptide_progs/**,

Name       | Description
-----------|----------------------
0_utils.sh | Code snippets
1_pgwas.sh | Association analysis.
2_meta_analysis.sh | Meta-analysis.
3_merge.sh   | Signal detection/classification, forest, Q-Q, Manhattan, LocusZoom, mean-by-genotype/dosage plots.
3.1_merge.sh | signal extraction
3.2_merge.sh | collection
3.3_merge.sh | graphical representation
-----------|----------------------

3.1-3.3 are spun off. In particular, CO3 and ITIH2 are resumed after the 12hr threshold for SLURM.

Prerequistes for a Manhattan/peptide association plot by `0_utils.sh` are

- a call to `gz()` for a compressed DR-filtered data; this is to be followed by execution of the R script inside.
- Ensembl-VEP (also step 3.2 above) but `ceuadmin/ensembl-vep/104` does not work on icelake; a remedy is made with `ceuadmin/ensembl-vep/111-icelake` however the loftee plugin is likely to require additional work -- as documented this is feasible with `vep -i variants.vcf --plugin LoFtool,scores_file.txt`.
