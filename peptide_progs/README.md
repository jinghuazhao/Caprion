# Peptide analysis

Peptide-level analysis mirrors protein-level analysis as in **peptide_progs/**,

Name       | Description
-----------|----------------------
0_utils.sh | Code snippets
1_pgwas.sh | Association analysis
2_meta_analysis.sh | Meta-analysis
3 | Signal identification
3.1_extract.sh | Signal extraction
3.2_collect.sh | Signal collection/classification
3.3_plot.sh | Forest, Q-Q, Manhattan, LocusZoom, mean-by-genotype/dosage plots

In particular, CO3 and ITIH2 are resumed after the 12hr threshold for SLURM.

Prerequistes for a Manhattan/peptide association plot by `0_utils.sh` are

- a call to `gz()` for a compressed DR-filtered data.
- Ensembl-VEP (also step 3.2 above) but `ceuadmin/ensembl-vep/104` has to be `ceuadmin/ensembl-vep/111-icelake` so is feasible with `vep -i variants.vcf --plugin LoFtool,scores_file.txt`.

---

The CSD3 icelake module `ceuadmin/R/4.3.3-icelake` now works as smoothly as `ceuadmin/R` for cclake
