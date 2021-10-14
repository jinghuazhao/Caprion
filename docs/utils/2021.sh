#!/usr/bin/bash

export dir=~/rds/projects/Caprion_proteomics
export caprion=${dir}/pilot
if [ ! -d ${caprion}/data3 ]; then mkdir ${caprion}/data3; fi
if [ ! -d ${caprion}/bgen3 ]; then mkdir ${caprion}/bgen3; fi

R --no-save -q < ${caprion}/utils/eSet.R

R --no-save -q < ${caprion}/utils/2021.R

R --no-save -q < ${caprion}/utils/UDP.R

pandoc 2021.md -o overlap.docx
