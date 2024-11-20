#!/usr/bin/bash

# Define the root path (assuming it was set earlier)
export root=${root:-/home/jhz22/Caprion/analysis/dup}

# Initialize directories
function initialise() {
  for d in sentinels/slurm json/gz means work qqmanhattanlz vep; do
    if [ ! -d ${root}/${d} ]; then mkdir -p ${root}/${d}; fi
  done
}

# Step 1: Create the task list
function step1_pqtl_list() {
  cat <<'EOL' > ${root}/dup-step1.sh
#!/usr/bin/bash

export root=${root}
export isotope=${1}  # Take isotope as an argument

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/R
module load samtools/1.13/gcc/zwxn7ug3
module load libiconv/1.16/intel/64iicvbf

export TMPDIR=${HPC_WORK}/work

function pgz() {
  # 1. extract all significant SNPs
  cat <(zcat ${root}/ZWK-1-${isotope}.fastGWA.gz | awk 'NR>1 && $7>=0.001 && $7 <= 0.999 && $10<=5e-8') \
      <(zcat ${root}/ZWK-1-${isotope}-chrX.fastGWA.gz | awk 'NR>1 && $7>=0.001 && $7 <= 0.999 && $10<=5e-8') | \
  sort -k1,1n -k3,3n | \
  gzip -f > ${root}/sentinels/${isotope}.p.gz
}

function _HLA()
# 2. handling HLA
{
  (
    zcat ${root}/ZWK-1-${isotope}.fastGWA.gz | awk -vOFS="\t" 'NR==1{$1="Chrom";$2="Start" "\t" "End";$3="SNP";print}'
    zcat ${root}/sentinels/${isotope}.p.gz | \
    awk -vOFS="\t" '{$1="chr"$1; snp=$2; $2=$3"\t"$3;$3=snp;print}' | \
    awk '!($1 == "chr6" && $2 >= 25392021 && $2 < 33392022)'
    zcat ${root}/sentinels/${isotope}.p.gz | \
    awk -vOFS="\t" '{$1="chr"$1; snp=$2; $2=$3"\t"$3;$3=snp;print}' | \
    awk '$1 == "chr6" && $2 >= 25392021 && $2 < 33392022' | \
      sort -k11,11g | \
      awk 'NR==1'
  ) > ${root}/sentinels/${isotope}.nold
  export lines=$(wc -l ${root}/sentinels/${isotope}.nold | cut -d' ' -f1)
  if [ $lines -eq 1 ]; then
     echo removing ${isotope}.nold with $lines lines
     rm ${root}/sentinels/${isotope}.nold
  fi
}

function sentinels() {
  # 3. merge and filter significant SNPs
  if [ -f ${root}/sentinels/${isotope}.nold ]; then
  (
    mergeBed -i ${root}/sentinels/${isotope}.nold -d 1000000 -c 11 -o min | \
    awk -v OFS="\t" -v trait=${isotope} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "trait"
      print $0, trait
    }'
  ) > ${root}/sentinels/${isotope}.merged
  (
    cut -f1-4,11 ${root}/sentinels/${isotope}.nold | \
    bedtools intersect -a ${root}/sentinels/${isotope}.merged -b - -wa -wb | \
    awk '$4==$10' | \
    awk -v OFS="\t" -v trait=${isotope} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "trait", "MarkerName", "CHR", "POS", "SNP", "P_check"
      chr=gsub(/chr/,"",$6)
      print $1,$2,$3,$4,$5,$6":"$7,chr,$7,$9,$10
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
  fi
}

# Now we run the functions
for cmd in pgz _HLA sentinels; do $cmd; done
EOL

  export N=$(head -1 ${root}/ZWK.pheno | awk '{print NF-2}')
  sed -i "s|ROOT|${root}|;s|_N_|${N}|" ${root}/dup-step1.sh
  chmod +x ${root}/dup-step1.sh
}

# Main script execution
export analysis=~/Caprion/analysis
export root=${analysis}/dup
initialise
step1_pqtl_list

# Now create the parallel execution environment:

# Extract the list of isotopes (one for each task), and run in parallel
isotopes=$(head -1 ${root}/ZWK.pheno | awk '{for(i=3;i<=NF;i++) print $i}')

# Use GNU Parallel to run the script in parallel
echo "$isotopes" | parallel -j 10 --bar ${root}/dup-step1.sh
