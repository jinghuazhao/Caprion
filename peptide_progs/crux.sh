#!/usr/bin/bash

module load ceuadmin/crux/4.1
export spectra=szwk901104i19801xms1

# 1: A Protein Database

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# 2: Data

# raw data --> mzML, done elsewhere
msconvert ${spectra}.raw --mzML

# 3: Perform Peptide Identification

# Run Crux with a chosen search engine Tide:
crux tide-index --overwrite T uniprot_sprot.fasta tide-index
crux tide-search --overwrite T --output-dir tide-output ${spectra}.mzML tide-index

# 4: Results

# Percolator validation
crux percolator --overwrite T --output-dir percolator-output tide-output/tide-search.target.txt

Rscript -e '
  to <- read.delim("tide-output/tide-search.target.txt")
  subset(to,grepl("PROC_HUMAN",protein.id))
  tp <- read.delim("percolator-output/percolator.target.peptides.txt")
  subset(tp,grepl("PROC_HUMAN",protein.id))
  psms <- read.delim("percolator-output/percolator.target.psms.txt")
  subset(psms,grepl("PROC_HUMAN",protein.id))
'

# utilities

crux version
crux get-ms2-spectrum ${spectra}.ms2
