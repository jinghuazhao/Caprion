#!/usr/bin/bash

function initialise()
{
  for d in sentinels/slurm fp means work qqmanhattanlz METAL/vep
  do
    if [ ! -d ${root}/${d} ]; then mkdir -p ${root}/${d}; fi
  done
}

function step1_pqtl_list()
{
cat <<'EOL'> ${root}/${protein}-step1.sb
#!/usr/bin/bash

#SBATCH --job-name=_1-LABEL
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem

#SBATCH --export ALL
#SBATCH --array=1-_N_
#SBATCH --output=ROOT/sentinels/slurm/_step1_%A_%a.o
#SBATCH --error=ROOT/sentinels/slurm/_step1_%A_%a.e

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/R/4.3.3-icelake
module load samtools/1.13/gcc/zwxn7ug3
module load libiconv/1.16/intel/64iicvbf

export TMPDIR=${HPC_WORK}/work
export protein=PROTEIN
export isotope=$(head -1 ${root}/${protein}.pheno | awk -vn=${SLURM_ARRAY_TASK_ID} '{print $(n+2)}')

function pgz()
# 1. extract all significant SNPs
{
# zcat METAL/{}-1.tbl.gz | head -1
  zcat ${root}/METAL/${isotope}-1.tbl.gz | \
  awk 'NR>1 && $6>=0.01 && $12<=log(5e-8)/log(10)' | \
  sort -k1,1n -k2,2n | \
  gzip -f > ${root}/sentinels/${isotope}.p.gz
}

function _HLA()
# 2. handling HLA
{
  (
    zcat ${root}/METAL/${isotope}-1.tbl.gz | awk -vOFS="\t" 'NR==1{$1="Chrom";$2="Start" "\t" "End";print}'
    zcat ${root}/sentinels/${isotope}.p.gz | \
    awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
    awk '!($1 == "chr6" && $3 >= 25392021 && $3 < 33392022)'
    zcat ${root}/sentinels/${isotope}.p.gz | \
    awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
    awk '$1 == "chr6" && $3 >= 25392021 && $3 < 33392022' | \
      sort -k13,13g | \
      awk 'NR==1'
  ) > ${root}/sentinels/${isotope}_nold.p
  export lines=$(wc -l ${root}/sentinels/${isotope}_nold.p | cut -d' ' -f1)
  if [ $lines -eq 1 ]; then
     echo removing ${isotope}_nold with $lines lines
     rm ${root}/sentinels/${isotope}_nold.p
  fi
}

function sentinels()
{
  (
    mergeBed -i ${root}/sentinels/${isotope}_nold.p -d 1000000 -c 13 -o min | \
    awk -v OFS="\t" -v trait=${isotope} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "trait"
      print $0, trait
    }'
  ) > ${root}/sentinels/${isotope}.merged
  (
    cut -f1-4,13 ${root}/sentinels/${isotope}_nold.p| \
    bedtools intersect -a ${root}/sentinels/${isotope}.merged -b - -wa -wb | \
    awk '$4==$10' | \
    cut -f1-5,9,10 | \
    awk -v OFS="\t" -v trait=${isotope} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "trait", "MarkerName", "CHR", "POS", "SNP", "P_check"
      chr=gsub(/chr/,"",$1)
      print $1,$2,$3,$4,trait,$1":"$3,chr,$3,$6,$7
    }'
  ) | uniq > ${root}/sentinels/${isotope}.sentinels

  Rscript -e '
    d <- file.path(Sys.getenv("root"),"sentinels")
    isotope <- Sys.getenv("isotope")
    f <- file.path(d,paste0(isotope,".sentinels"))
    m <- read.table(f,header=TRUE,as.is=TRUE)
    dim(m)
    head(m)
    suppressMessages(library(dplyr))
    t <- m %>% group_by(trait,Chrom,Start,End) %>% slice(which.min(P))
    t
    P <- with(m,P)
    p <- table(P)[table(P)>1]
    print(p)
    m <- subset(t,MarkerName!=".")
    cols <- c(1:5,9)
    signal_file <- file.path(d,paste0(isotope,".signals"))
    if (exists(signal_file)) unlink(signal_file)
    write.table(m[,cols],file=signal_file,row.names=FALSE,quote=FALSE,sep="\t")
  '
}

for cmd in pgz _HLA sentinels; do $cmd; done
EOL

export N=$(head -1 ${root}/${protein}.pheno | awk '{print NF-2}')
sed -i "s|ROOT|${root}|;s|LABEL|${protein}|;s|PROTEIN|${protein}|;s|_N_|${N}|" ${root}/${protein}-step1.sb
}

source 3.0_merge.sh

# all proteins:
xargs -n 2 < ${analysis}/peptide_progs/benchmark2.names | \
grep -n -f ${analysis}/peptide_progs/benchmark2.names -w ${varlist} | \
while IFS=":" read -r protein_index protein; do
    export protein_index
    export protein
    echo ${protein_index} ${protein}
    export root=~/Caprion/analysis/peptide/${protein}
    export pheno=${analysis}/peptide/${protein}/${protein}.pheno
    export N=$(awk 'NR==1{print NF-2}' ${pheno})
    initialise
  # 1. Extraction of signals
    echo Step 1.
    step1_pqtl_list
    sbatch ${root}/${protein}-step1.sb
done
