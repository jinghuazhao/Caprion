#!/usr/bin/bash
# https://crux.ms/fileformats.html

module load ceuadmin/wine/8.21
module load ceuadmin/pwiz/3_0_24163_9bfa69a-wine
export ZWK=~/Caprion/pre_qc_data/spectral_library_ZWK
export raw=szwk901104i19801xms1

# option 1. singularity
singularity --version
export SIF=pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif
if [ ! -f "${f}" ] && [ -x "${f}" ]; then
   singularity pull --name ${SIF} docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
 # http://localhost:8080, Admin, password
   ln -sf /rds/project/rds-4o5vpvAowP0/software/.apptainer/ ${HOME}/.apptainer
   singularity pull docker://quay.io/galaxy/introduction-training
fi
for format in --mzML --mzXML --mz5 --mzMLb --mgf --text --ms1 --cms1 --ms2 --cms2
do
singularity exec --env WINEDEBUG=-all \
                  -B /rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics/spectral_library_ZWK/:/data \
                      ${SIF} \
                      wine msconvert ${format} /data/${raw}.raw
done

# option 2. ThermoRawFileParser (ThermoFileReader)
## CLI
module load ceuadmin/ThermoRawFileParser
cd $ThermoRawFileParser_HOME;
ThermoRawFileParser.exe -i $ZWK/$raw
cd -

## GUI
module load ceuadmin/ThermoRawFileParserGUI
cd $ThermoRawFileParser_HOME/resources/ThermoRawFileParser
java -jar $ThermoRawFileParser_HOME/ThermoRawFileParser-1.7.4.jar
cd -

# option 3 (legacy)
wine64 $(which msconvert.exe) --mgf ${ZWK}/${raw}
wine64 $(which msconvert.exe) --mzML ${ZWK}/${raw}

export sprot=~/rds/public_databases/UniProt/uniprot_sprot.fasta.gz
export sprot=uniprot_sprot.fasta
export uniref100=~/rds/public_databases/UniProt/uniref100.fasta.gz
export uniprot='Human_database_including_decoys_(cRAP_added)'

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

# paket add Mono.Unix --version 7.1.0-final.1.21458.1
module load mono/5.0.1.1
module load ceuadmin/dotnet/6.0.423
module load ceuadmin/icu/50.2 ceuadmin/OpenMS/3.0.0-pre-develop-2022-09-28
module load ceuadmin/tandem/2017.2.1.4

# decoy database
DecoyDatabase -in ${uniprot}.fasta -out decoy_${uniprot}.fasta

# Convert raw data to mzML format using ThermoRawFileParser.exe
# msconvert qExactive01819.raw qExactive01819.mzML
export spectra=szwk901104i19801xms1
ln -sf /rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics/spectral_library_ZWK/${spectra}.raw rawdata.raw
FileConverter -in rawdata.raw -out rawdata.mzML

## 1.
# peptide identification & more showcase of singularity
XTandemAdapter -in ${spectra}.mzML -database ${unprot}.fasta -xtandem_executable $(which tandem.exe)

# Perform peak picking
PeakPickerHiRes -in ${spectra}.mzML -out ${spectra}_picked.mzML

# Run peptide identification using Comet
CometAdapter -in ${spectra}_picked.mzML -out ${spectra}_comet.idXML -database ${uniprot}.fasta

# Convert identification results to OpenMS format using comet.exe
IDFileConverter -in ${spectra}_comet.idXML -out comet.idXML

# Map peptide identifications to features
IDMapper -in ${spectra}_picked.mzML -id comet.idXML -out mapped.mzML

# Annotate peptides with protein information
PeptideIndexer -in comet.idXML -fasta ${uniprot}.fasta -out indexed.idXML

singularity exec -B /usr/local/Cluster-Apps/ceuadmin/OpenMS/3.0.0-pre-develop-2022-09-28/bin/:/openms \
                 -B /usr/local/Cluster-Apps/ceuadmin/tandem/2017.2.1.4/:/bin \
                 -B /rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/OpenMS/tutorials:/data \
            introduction-training_latest.sif /openms/XTandemAdapter \
                -in /data/qExactive01819_profile.mzml -database /data/'Human_database_including_decoys_(cRAP_added).fasta'
## 2.

export spectra2=szwk021704i19101xms1
export db=Human_database_including_decoys_(cRAP_added).fasta

# Feature detection: 16193 features found.
FeatureFinderCentroided -in ${spectra2}.mzML -out ${spectra2}.featureXML

# Peptide Identification
XTandemAdapter -in ${spectra2}.mzML -out identifications.idXML -database ${db}

# Validation and Filtering
IDFilter -in identifications.idXML -out filtered_identifications.idXML -score:pep 0.01

# Peptide Indexing
PeptideIndexer -in filtered_identifications.idXML -out indexed_identifications.idXML -fasta ${db}

# Protein Inference
ProteinQuantifier -in indexed_identifications.idXML -out protein_quant.csv -consensus:protein_level

## 3.

INPUT_DIR="."
OUTPUT_DIR="."

for file in $INPUT_DIR/sz*mzML
do
    filename=$(basename "$file" .mzML)

  # Perform peak detection
    PeakPickerHiRes -in "$file" -out "$OUTPUT_DIR/${filename}_peaks.mzML"

  # Detect and quantify features
    FeatureFinderCentroided -in "$OUTPUT_DIR/${filename}_peaks.mzML" -out "$OUTPUT_DIR/${filename}.featureXML"
done

# Align features across samples (if applicable)
FeatureLinkerUnlabeledQT -in $OUTPUT_DIR/*featureXML -out $OUTPUT_DIR/aligned.consesusXML

# Statistical analysis
ProteinQuantifier -in aligned.consensusXML -out quantified.csv
}

function docker()
# under Ubuntu 22.04 LTS /home/$USER/D/Downloads:/data
{
  docker run -it --rm -e WINEDEBUG=-all -v ${ZWK}:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses \
         wine msconvert /data/szwk901104i19801xms1.raw --filter "peakPicking true 1-"
  docker run -p 8080:80 quay.io/galaxy/introduction-training
  module load ceuadmin/docker/24.0.5
}
