#!/usr/bin/bash

export rt="caprion"
cut -f6 ${rt}.merge | sed '1d' | sort | uniq > ${rt}.snp

function ps()
{
  module load ceuadmin/phenoscanner
  split -l 500 --numeric-suffix ${rt}.snp ${rt}.snp.
  for i in `seq 0 4`
  do
    echo ${rt}.snp.0${i}.pQTL
    phenoscanner -s T -c pQTL -x EUR -p 0.0000001 --r2 0.6 -i ${rt}.snp.0${i} -o ${rt}.snp.0${i}.pQTL
  done
  (
    for i in `seq 0 4`
    do
      if [ $i -eq 0 ]; then
         cat ${rt}.snp.0${i}.pQTL
      else
         sed '1d' ${rt}.snp.0${i}.pQTL
      fi
      rm ${rt}.snp.0${i}.pQTL
    done
  ) > ${rt}.pQTL
}
