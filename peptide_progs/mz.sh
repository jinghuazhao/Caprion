#!/usr/bin/bash
# https://crux.ms/fileformats.html

module load ceuadmin/wine/8.21
module load ceuadmin/pwiz/3_0_24163_9bfa69a-wine
export ZWK=~/Caprion/pre_qc_data/spectral_library_ZWK
export raw=szwk901104i19801xms1

singularity --version
export SIF=pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif
if [ ! -f "${f}" ] && [ -x "${f}" ]; then
   singularity pull --name ${SIF} docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
fi
for format in --mzML --mzXML --mz5 --mzMLb --mgf --text --ms1 --cms1 --ms2 --cms2
do
singularity exec --env WINEDEBUG=-all \
                  -B /rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics/spectral_library_ZWK/:/data \
                      ${SIF} \
                      wine msconvert ${format} /data/${raw}.raw
done

module load ceuadmin/crux/4.1
export sprot=~/rds/public_databases/UniProt/uniprot_sprot.fasta.gz
export sprot=uniprot_sprot.fasta
export uniref100=~/rds/public_databases/UniProt/uniref100.fasta.gz
export mgf=${raw}.mgf
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

function boxcar()
{
  module load ceuadmin/Anaconda3/2023.09-0
  pip install pyteomics
  python BoxCar.py
  python pyteomics.py
}

function openms()
{
# A Workflow for Peptide Identification

##   Preprocessing Raw Data:
##       Convert raw MS data to an appropriate format using FileConverter.
##       Perform peak picking and feature detection using FeatureFinderCentroided.

##   Peptide Identification:
##      Run peptide identification using search engines like Comet, XTandem, or MSGFPlus.
##      Convert the identification results to OpenMS format using IDFileConverter.

##   Mapping and Annotation:
##      Map the identified peptides to detected features using IDMapper.
##      Annotate peptide sequences with protein information using PeptideIndexer.

##  Quantification and Analysis:
##      Quantify peptides and proteins using ProteinQuantifier or other relevant tools.
##      Perform further statistical analysis and visualization using additional OpenMS tools or external software.

module load ceuadmin/icu/50.2 ceuadmin/OpenMS/3.0.0-pre-develop-2022-09-28
module load ceuadmin/tandem/2017.2.1.4
export spectra=szwk901104i19801xms1
ln -sf /rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics/spectral_library_ZWK/${spectra}.raw rawdata.raw

# Convert raw data to mzML format using ThermoRawFileParser.exe
FileConverter -in rawdata.raw -out rawdata.mzML

# peptide identification
XTandemAdapter -in ${spectra}.mzML -database uniprot_sprot.fasta -xtandem_executable $(which tandem.exe)

# Perform peak picking
PeakPickerHiRes -in ${spectra}.mzML -out picked.mzML

# Run peptide identification using Comet
CometAdapter -in picked.mzML -out comet.idXML -database uniprot_sprot.fasta

# Convert identification results to OpenMS format using comet.exe
IDFileConverter -in comet.idXML -out comet.idXML

# Map peptide identifications to features
IDMapper -in picked.mzML -id comet.idXML -out mapped.mzML

# Annotate peptides with protein information
PeptideIndexer -in comet.idXML -fasta uniprot.fasta -out indexed.idXML
}
