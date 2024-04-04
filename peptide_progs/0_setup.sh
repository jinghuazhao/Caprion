export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export signals=${analysis}/work/caprion${suffix}.signals
export varlist=${analysis}/output/caprion${suffix}.varlist

function with_pQTL_only1()
# only those with pQTLs:
{
  export n_with_signals=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | wc -l)
  for i in $(echo $(seq ${n_with_signals} | grep -w -f <(sed 's/, /\n/g' benchmark2.lst)))
  do
    export signal_index=${i}
    export signal_list=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | awk 'NR==ENVIRON["signal_index"]' | grep -f - -n -w ${signals} | cut -d':' -f1)
    export protein_index=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | awk 'NR==ENVIRON["signal_index"]' | grep -f - -n -w ${varlist} | cut -d':' -f1)
    export protein=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | awk 'NR==ENVIRON["signal_index"]')
    echo $protein, $signal_index, $protein_index, $signal_list
    sb ${protein_index}
  done
}
# 79 proteins have errors since they took longer than 12hrs to finish:
# for i in $(grep error ${analysis}/peptide/*/*.e | sed 's|/|\t|g' | cut -f7 | grep -n -w -f - ${pilot}/work/caprion.varlist | cut -d':' -f1)
# Some batches contain no data
# 72 160 237 382 for BROX, CT027, GHRL, NCF2

function with_pQTL_only2()
# only those with pQTLs
{
  export n_with_signals=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | wc -l)
  for i in $(echo $(seq ${n_with_signals} | grep -w -f <(sed 's/, /\n/g' benchmark2.lst)))
  do
    export signal_index=${i}
    export protein=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | awk 'NR==ENVIRON["signal_index"]')
    export root=~/Caprion/analysis/peptide/${protein}
    export pheno=${analysis}/peptide/${protein}/${protein}.pheno
    export N=$(awk 'NR==1{print NF-2}' ${pheno})
    echo ${signal_index}, ${protein}
    if [ ! -d ${root}/METAL ]; then mkdir ${root}/METAL; fi
    METAL_list
    METAL_files
    METAL_analysis_sbatch
  done
}

function with_pQTL_only3()
# only those with pQTLs
{
export n_with_signals=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | wc -l)
for i in $(echo $(seq ${n_with_signals} | grep -w -f <(sed 's/, /\n/g' benchmark.lst)))
do
  export signal_index=${i}
  export protein=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | awk 'NR==ENVIRON["signal_index"]')
  export root=${analysis}/peptide/${protein}
  export pheno=${root}/${protein}.pheno
  export N=$(awk 'NR==1{print NF-2}' ${pheno})
  export all_peptides=$(head -1 ${pheno} | cut -d' ' -f1,2 --complement)
  export dir=${root}/qqmanhattanlz
  echo ${signal_index}, ${protein}
  initialise
# 1. Extraction of signals
  echo Step 1:
  step1_pqtl_list
  sbatch ${root}/${protein}-step1.sb
# 2. Collection of signals
  echo Step 2:
  step2_pqtl_collect
# fplz should be here
  sbatch ${root}/${protein}-step2.sb
# 3. Graphical representation
  echo Step 3:
  export pqtl_peptides=$(sed '1d' ${root}/${protein}.signals | cut -f1 | sort -k1,1n | uniq)
  export array=$(grep -n -f <(echo ${pqtl_peptides} | tr ' ' '\n') <(echo ${all_peptides} | tr ' ' '\n') | cut -d':' -f1 | tr '\n' ',' | sed 's/.$//')
# 12hr-timeout proteins
# export i=129, 299 for CO3, ITIH2; 47, 179, 184, 492 for APOB, EPCR, ERAP2, PROC
# export todo_peptides=$(ls ${analysis}/peptide/CO3/fp/*pdf | xargs -l basename -s -fp.pdf | sort -k1,1 | uniq | grep -f - -v <(sed 's/ /\n/g' <(echo ${pqtl_peptides})))
# export array=$(grep -n -f <(echo ${todo_peptides} | tr ' ' '\n') <(echo ${all_peptides} | tr ' ' '\n') | cut -d':' -f1 | tr '\n' ',' | sed 's/.$//')
  for batch in {1..3}
  do
  (
    head -1 ${root}/${protein}.pheno
    grep -f <(cut -d' ' -f1 ${pilot}/output/caprion-${batch}.id) ${root}/${protein}.pheno
  ) > ${root}/work/${protein}-${batch}.pheno
  done
  step3_pqtl_summary
  sbatch ${root}/${protein}-step3.sb
# pdfjam ${dir}/*_qqmanhattan.pdf --nup 1x1 --landscape --papersize '{7in,12in}' --outfile ${root}/qq+manhattan.pdf
# qpdf --empty --pages $(ls ${dir}/*_qqmanhattan.pdf) -- ${root}/qq+manhattan.pdf
done
}
