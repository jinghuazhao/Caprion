#!/usr/bin/bash

module load ceuadmin/crux/4.1
module load ceuadmin/pwiz/3_0_24163_9bfa69a-wine
module load ceuadmin/wine/8.21
export sprot=~/rds/public_databases/UniProt/uniprot_sprot.fasta.gz
export sprot=uniprot_sprot.fasta
export uniref100=~/rds/public_databases/UniProt/uniref100.fasta.gz

export raw=szwk901104i19901xms1.raw
# https://crux.ms/fileformats.html
# module load ceuadmin/pwiz/3_0_24156_80747de not fully functional
# only works on Windows
{
  export ZWK=~/Caprion/pre_qc_data/spectral_library_ZWK
  export mgf=$(echo ${raw} | sed 's/raw/mgf/')
  export out=$(echo ${raw} | sed 's/raw/txt/')
  echo ${mgf}
  if [ ! -f ${mgf} ]; then
     wine64 $(which msconvert.exe) --mgf ${ZWK}/${raw}
  fi
}

R --no-save < crux.R

function yeast()
# Rscript multidecoy_peptides.R yeast_parameters.txt
{
  wget https://bitbucket.org/noblelab/multi-competition-fdr/raw/195d3608e427b17eebe555ca280c6db42b52678c/multidecoy_peptides.R
  wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/04/PXD002726/KFSwYeastCrTry_150506.fasta
  wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/04/PXD002726/Yeast_In-gel_digest_2.mgf
}
