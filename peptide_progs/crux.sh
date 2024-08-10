#!/usr/bin/bash

module load ceuadmin/crux/4.1
export spectra=szwk021704i19101xms1

# 1: A Protein Database

# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
# gunzip uniprot_sprot.fasta.gz

# 2: Data

# raw data --> mzML, done elsewhere
# msconvert ${spectra}.raw --mzML

# 3: Perform Peptide Identification

# Run Crux with a chosen search engine Tide:
crux tide-index --overwrite T uniprot.fasta tide-index
crux tide-search --overwrite T --output-dir tide-output ${spectra}.mzML tide-index

# 4: Results

# Percolator validation
crux percolator --overwrite T --output-dir percolator-output tide-output/tide-search.target.txt

Rscript -e '
  to <- read.delim("tide-output/tide-search.target.txt")
  tp <- read.delim("percolator-output/percolator.target.peptides.txt")
  psms <- read.delim("percolator-output/percolator.target.psms.txt")
  sink("PROC-crux.md")
  knitr::kable(subset(to,grepl("PROC_HUMAN",protein.id)),caption="tide-search.target")
  knitr::kable(subset(tp,grepl("PROC_HUMAN",protein.id)),caption="percolator.target.peptides")
  knitr::kable(subset(psms,grepl("PROC_HUMAN",protein.id)),caption="percolator.target.psms")
  sink()
'

# utilities

crux version
crux get-ms2-spectrum ${spectra}.ms2
