#!/usr/bin/bash

function pgz1()
# 1. extract all significant SNPs
{
  ls bgen/*.gz | grep invn | \
  sed 's|bgen/||g;s/-plink2//g;s/.gz//g' | \
  parallel -j6 -C' ' '
  (
    zcat bgen/{}-plink2.gz | head -1
    zcat bgen/{}-plink2.gz | awk "NR>1 && \$12 <=1e-5" | sort -k1,1n -k2,2n
  ) | gzip -f > sentinels/{}.p.gz'
}

function _HLA1()
# 2. handling HLA
{
  for p in $(ls bgen/*.gz | grep invn | sed 's|bgen/||g;s/-plink2//g;s/.gz//g')
  do
    (
      zcat bgen/${p}-plink2.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      (
        zcat sentinels/${p}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk 'NR>1 && !($1 == "6" && $3 >= 25392021 && $3 < 33392022)'
        zcat sentinels/${p}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk '$1 == "6" && $3 >= 25392021 && $3 < 33392022' | \
        sort -k13,13g | \
        awk 'NR==1'
      ) | \
      sort -k1,1n -k2,2n -k3,3n | \
      awk -v OFS="\t" '{$1="chr" $1};1'
    ) > sentinels/${p}${tag}.p
    export lines=$(wc -l sentinels/${p}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}${tag} with $lines lines
      rm sentinels/${p}${tag}.p
    fi
  done
}

function pgz2()
# 1. extract all significant SNPs
{
  ls ${src}/*.gz | grep invn | \
  sed "s|${src}/||g;s/-plink2//g;s/.gz//g" | \
  parallel -j6 -C' ' --env src --env dest --env p '
  (
    zcat ${src}/{}-plink2.gz | head -1
    zcat ${src}/{}-plink2.gz | awk -v p=${p} "NR>1 && \$12 <=p" | sort -k1,1n -k2,2n
  ) | gzip -f > ${dest}/{}.p.gz'
}

function _HLA2()
# 2. handling HLA
{
  for p in $(ls ${src}/*.gz | grep invn | sed "s|${src}/||g;s/-plink2//g;s/.gz//g")
  do
    (
      zcat ${src}/${p}-plink2.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      (
        zcat ${dest}/${p}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk 'NR>1 && !($1 == "6" && $3 >= 25392021 && $3 < 33392022)'
        zcat ${dest}/${p}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk '$1 == "6" && $3 >= 25392021 && $3 < 33392022' | \
        sort -k13,13g | \
        awk 'NR==1'
      ) | \
      sort -k1,1n -k2,2n -k3,3n | \
      awk -v OFS="\t" '{$1="chr" $1};1'
    ) > ${dest}/${p}${tag}.p
    export lines=$(wc -l ${dest}/${p}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}${tag} with $lines lines
      rm ${dest}/${p}${tag}.p
    fi
  done
}

function pgz2group()
# 1. extract all significant SNPs
{
  ls ${src}/*.gz | grep invn | grep ${group} | \
  sed "s|${src}/||g;s/-${group}-plink2//g;s/.gz//g" | \
  parallel -j6 -C' ' --env src --env dest --env p --env group '
  (
    zcat ${src}/{}-${group}-plink2.gz | head -1
    zcat ${src}/{}-${group}-plink2.gz | awk -v p=${p} "NR>1 && \$12 <=p" | sort -k1,1n -k2,2n
  ) | gzip -f > ${dest}/{}-${group}.p.gz'
}

function _HLA2group()
# 2. handling HLA
{
  for p in $(ls ${src}/*.gz | grep invn | grep ${group} | sed "s|${src}/||g;s/-${group}-plink2//g;s/.gz//g")
  do
    (
      zcat ${src}/${p}-${group}-plink2.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      (
        zcat ${dest}/${p}-${group}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk 'NR>1 && !($1 == "6" && $3 >= 25392021 && $3 < 33392022)'
        zcat ${dest}/${p}-${group}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk '$1 == "6" && $3 >= 25392021 && $3 < 33392022' | \
        sort -k13,13g | \
        awk 'NR==1'
      ) | \
      sort -k1,1n -k2,2n -k3,3n | \
      awk -v OFS="\t" '{$1="chr" $1};1'
    ) > ${dest}/${p}-${group}${tag}.p
    export lines=$(wc -l ${dest}/${p}-${group}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}-${group}${tag} with $lines lines
      rm ${dest}/${p}-${group}${tag}.p
    fi
  done
}

export tag=_nold
export src=bgen2
export dest=sentinels
export p=5e-8
if [ ! -d ${dest} ]; then mkdir -p ${dest}; fi
# for group in group1 group2; do export group=${group}; pgz2group; _HLA2group; done

pgz2
_HLA2
