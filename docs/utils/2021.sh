#!/usr/bin/bash

export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot
if [ ! -d ${caprion}/data3 ]; then mkdir ${caprion}/data3; fi
if [ ! -d ${caprion}/bgen3 ]; then mkdir ${caprion}/bgen3; fi

function symlink()
{
  cd ${caprion}
  for f in INTERVALdata_17FEB2021.csv INTERVAL_OmicsMap_20210217.csv pilotsMap_17FEB2021.csv; do ln -sf ../$f; done
  for f in interval_caprion_pilot_samples_phenotype_data.tsv ZWK_EDR_20191002.xlsx ZWK_Summary\ Report_20191018.pptx; do ln -sf ../$f ; done
}

function run_analysis()
{
  R --no-save -q < ${caprion}/utils/eSet.R
  R --no-save -q < ${caprion}/utils/2021.R
  R --no-save -q < ${caprion}/utils/UDP.R
  sbatch ${caprion}/utils/qctool.sb
  ${caprion}/utils/qctool.sh
  sbatch ${caprion}/utils/plink2.sb
}

function overlap()
{
  export threshold=$1
  export a=${caprion}/bgen2/${threshold}/caprion-invn.sentinels
  export b=${caprion}/bgen3/protein/${threshold}/caprion-invn.sentinels
  bedtools intersect -a <(awk '{if(NR>1) $1="chr"$1;print}' $a | tr ' ' '\t') \
                     -b <(awk '{if(NR>1) $1="chr"$1;print}' $b | tr ' ' '\t') \
                     -wa -wb -loj | \
  awk '$8!="." && $5==$12' | \
  wc -l

  bedtools intersect -a <(sed '1d;s/_All_invn//g' $a | \
                          awk '{$1="chr"$1;print}' | \
                          tr ' ' '\t') \
                     -b <(sed '1d;s/protein-//;s/_invn//' $b | \
                          awk '{$1="chr"$1;print}' | \
                          tr ' ' '\t') \
                     -wa -wb -loj | \
  awk '$8!="." && $5==$12'
}

function run_overlap()
{
  for p in 1e-5 5e-8
  do
    echo ${p}
    echo All overlaps
    overlap ${p} | wc -l
    echo Replicated pQTLs
    overlap ${p} | cut -f1-7 | uniq | wc -l
  done
  overlap ${p} | cut -f2-3,9-10 --complement | sed 's/chr//g;s/[0-9]*://g' | \
  awk '{$1=$1":"$4;$6=$6":"$9};1' | cut -d' ' -f4,9 --complement | \
  awk '{print $1,$2,$3"-"$4,$5,$6,$7"-"$8}' | tr ' ' '|'
  overlap ${p} | awk '(NR>1&&$7==$14){print $5,$NF}' | \
  parallel -C' ' --env caprion '
    zgrep -w {2} ${caprion}/bgen2/{1}_All_invn-plink2.gz | \
    awk -v protein={1} -v OFS="|" "{print protein\"-\"\$3,\$4\"/\"\$5,\$9,\$10,\$12}"
    zgrep -w {2} ${caprion}/bgen3/protein-{1}_invn-plink2.gz | \
    awk -v protein={1} -v OFS="|" "{print protein\"-\"\$3,\$4\"/\"\$5,\$9,\$10,\$12}"
  '
}

pandoc 2021.md -o 2021.docx
