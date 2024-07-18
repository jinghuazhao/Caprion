#!/usr/bin/bash
# https://crux.ms/fileformats.html

module load ceuadmin/wine/8.21
module load ceuadmin/pwiz/3_0_24163_9bfa69a-wine
export ZWK=~/Caprion/pre_qc_data/spectral_library_ZWK
export raw=szwk901104i19801xms1.raw

singularity --version
export SIF=pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif
if [ ! -f "${f}" ] && [ -x "${f}" ]; then
   singularity pull --name ${SIF} docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
fi
singularity exec --env WINEDEBUG=-all \
                  -B /rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics/spectral_library_ZWK/:/data \
                      ${SIF} \
                      wine msconvert /data/${raw} --filter "peakPicking true 1-"

module load ceuadmin/crux/4.1
export sprot=~/rds/public_databases/UniProt/uniprot_sprot.fasta.gz
export sprot=uniprot_sprot.fasta
export uniref100=~/rds/public_databases/UniProt/uniref100.fasta.gz
R --no-save < crux.R

function legacy()
{
  export mgf=$(echo ${raw} | sed 's/raw/mgf/')
  export mzML=$(echo ${raw} | sed 's/raw/mzML/')
  if [ ! -f ${mgf} ]; then
     echo ${mgf}
     wine64 $(which msconvert.exe) --mgf ${ZWK}/${raw}
  fi
  if [ ! -f ${mzML} ]; then
     echo ${mzML}
     wine64 $(which msconvert.exe) --mzML ${ZWK}/${raw}
  fi
}

function docker()
# under Ubuntu 22.04 LTS /home/$USER/D/Downloads:/data
{
  docker run -it --rm -e WINEDEBUG=-all -v ${ZWK}:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses \
         wine msconvert /data/szwk901104i19801xms1.raw --filter "peakPicking true 1-"
  module load ceuadmin/docker/24.0.5
}

function yeast_files()
# Rscript multidecoy_peptides.R yeast_parameters.txt
{
  wget https://bitbucket.org/noblelab/multi-competition-fdr/raw/195d3608e427b17eebe555ca280c6db42b52678c/multidecoy_peptides.R
  wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/04/PXD002726/KFSwYeastCrTry_150506.fasta
  wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/04/PXD002726/Yeast_In-gel_digest_2.mgf
}
