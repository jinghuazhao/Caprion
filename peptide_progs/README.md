# Peptide analysis

Peptide-level analysis

This mirrors protein-level analysis as in peptide_progs/,

```
0_utils.sh  
1_pgwas.sh
2_meta_analysis.sh
3_merge.sh
```

with respect to 

1. Association analysis.
2. Meta-analysis.
3. Signal detection/classification, forest, Q-Q, Manhattan, LocusZoom, mean-by-genotype/dosage plots.

NB. `3_merge.sh` itself is comprised of three steps on signal extraction, collection, and graphical representation which is operationally done through commenting/uncommenting calls to each step exclusively. In particular, CO3 and ITIH2 are resumed after the 12hr threshold for SLURM.
