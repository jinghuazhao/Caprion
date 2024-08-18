#!/bin/bash

function sb()
{
cat << 'EOL' > ${sbatch}
#!/bin/bash

#SBATCH --job-name=_crux-search
#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4  # Use 4 cores per task

#SBATCH --output=ANALYSIS/crux/RAW.o
#SBATCH --error=ANALYSIS/crux/RAW.e

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/crux/4.1
export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=ANALYSIS
export PERL5LIB=
export raw=RAW

function crux_search()
# database search with Tide:
{
#1. creating a peptide index file from human proteins
  crux tide-index \
       --overwrite T --max-mods 3 --missed-cleavages 3 FASTA human-idx

#2. searching a spectrum dataset (UPS1.mzML.gz) with tide-search, large RAM required
  crux tide-search \
       --overwrite T --concat T --top-match 1 --num-threads 4 --precursor-window 100 --mz-bin-width 1.0005079 \
       --precursor-window-type ppm --use-tailor-calibration T \
       --output-dir results ${raw}.mzML.gz human-idx

#3. running Percolator:
  crux percolator --overwrite T --output-dir results results/tide-search.txt
}

cd ${analysis}/crux
crux_search
EOL
sed -i "s|ANALYSIS|${analysis}|;s|FASTA|${fasta}|;s|RAW|${raw}|" ${sbatch}
sbatch ${sbatch}
}

export analysis=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis
export raw=szwk901104i19801xms1
export raw=12March2019-Lumos-DIA-HuAD-Hipp-B1-ind-8mz-ovlp-400to1000-HZR03
export fasta=uniprot-proteome_UP000005640+reviewed_yes.fasta
export fasta=uniprot.fasta
export sbatch=${analysis}/crux/commands/${raw}.sb
sb

function test()
{
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
  crux assign-confidence --overwrite T --estimation-method tdc --output-dir assign-confidence \
                         --decoy-prefix decoy_ percolator-output/percolator.target.psms.txt
# utilities

  crux version
  crux generate-peptides --enzyme trypsin --missed-cleavages 2 --decoy-format peptide-reverse uniprot.fasta
  crux get-ms2-spectrum ${spectra}.ms2

  export uniprot=uniprot
  export mgf=${raw}
  module load ceuadmin/crux/4.1
  R --no-save < crux.R
}
