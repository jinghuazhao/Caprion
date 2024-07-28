!/bin/bash

cat << 'EOL' > ${sbatch}
#!/bin/bash

#SBATCH --job-name=_ZWK
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
#1. creating a peptide index file from human proteins (uniprot-proteome_UP000005640+reviewed_yes.fasta),
  crux tide-index \
       --overwrite T --max-mods 3 --mods-spec 2M+15.9949,2STY+79.966331 --missed-cleavages 3 \
       uniprot-proteome_UP000005640+reviewed_yes.fasta human-idx

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

export raw=szwk901104i19801xms1
export analysis=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis
export sbatch=${analysis}/crux/${RAW}.sbatch
sed -i "s/ANALYSIS/${analysis};s/RAW/${raw}/" ${sbatch}

sbatch ${sbatch}
