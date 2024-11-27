# Miscellaneous analysis

This section accommodates many largely independent tasks.

Implementation might well be generic so that both proteins and peptides are covered.

## Programs and applications

These are summarised in the following table,

Program   | Description
----------|------------------------------------------------------------------------
`csq.sh` | Consequences of variants
`Caprion_deCODE_UKB_PPP.sh` | Caprion/deCODE/UKB-PPP replication
`eSet.sh` | ExpresssionSet implementations
`glmnet_pense.sh` | glmnet/pense modeling
`impute.sb` | imputation experiments
`json.sh` | JSON file generation
`peptideAssociationPlot.sh` | protein Manhattan-peptide signal plots
`dup-pgwas.sh` | pGWAS for duplicated proteins
`dup-extract.sh` | pQTL extractions
`dup-json.sh` | LocusZoom.js plots
`dup-plot.sh` | pQTL plots
`pqtlGWAS.R` | pQTL-GWAS lookup
`tables.sh` | `Supplementary-Tables.xlsx` generator
`ToDo.sh` | various staged experiments

NB: `impute.sb` employs `impute_parallel()` when N(isotope groups) > 500. Nevertheless,
when coming to protein requantification this is an option to use the orginal intensity
data.

## Legacy codes

- `compare.sb`. earlier contrast with deCODE/UKB-PPP.
- `inf1.sh`. snapshot from SCALLOP/INF meta-analysis.

Created on **20/11/2024**
