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

R --no-save -q < ${caprion}/utils/eSet.R
R --no-save -q < ${caprion}/utils/2021.R
R --no-save -q < ${caprion}/utils/UDP.R

pandoc 2021.md -o overlap.docx

# --- legacy ---

sed '1d' ${caprion}/data3/UDP.tsv | grep -n NA | cut -f1 | sed 's/:NA//' | tr '\n' ' ' > ${caprion}/data3/na_UDP.lst

