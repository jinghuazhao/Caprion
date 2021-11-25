export caprion=~/Caprion/pilot
export SLURM_ARRAY_TASK_ID=1

function phase3()
{
  export data=${1}
  export col=$(head -1 ${caprion}/data3/${data}.pheno | sed 's/ /\n/g' | grep "_invn$" | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]')
  for v in ${col};
  do
      echo ${v}
      regenie --step 1 --force-step1 \
              --bgen ${caprion}/data3/caprion-1.bgen --sample ${caprion}/data3/${data}.sample \
              --phenoFile ${caprion}/data3/${data}.pheno --phenoCol ${v} \
              --qt --test additive --bsize 10000 \
              --out ${caprion}/work/${data}-${v}
      regenie --step 2 \
              --bgen ${caprion}/data3/caprion.01.bgen --sample ${caprion}/data3/${data}.sample \
              --phenoFile ${caprion}/data3/${data}.pheno --phenoCol ${v} \
              --qt --test additive --gz --bsize 20000 --pred ${caprion}/work/${data}-${v}.list \
              --out ${caprion}/work/${data}-${v}-begenie
  done
}

phase3 protein
