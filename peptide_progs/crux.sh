#!/usr/bin/bash

module load ceuadmin/crux/4.1
export spectra=szwk901104i19801xms1

# 1: A Protein Database

Download a FASTA file of a protein database, for example, from UniProt.

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# 2: Data

# raw data --> mzML
msconvert ${spectra}.raw --mzML

# 3: Perform Peptide Identification

# Run Crux with a chosen search engine Tide:
crux tide-index uniprot_sprot.fasta tide-index
crux tide-search --output-dir tide-output ${spectra}.mzML tide-index

# 4: Results

# Percolator validation
crux percolator --output-dir percolator-output tide-output/tide-search.target.txt
# Generate a report
crux q-ranker --output-dir qranker-output percolator-output/percolator.target.psms.txt
