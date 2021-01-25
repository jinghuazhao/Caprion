#!/usr/bin/bash

export caprion=${HOME}/Caprion
export protein=EPCR

qpdf ${caprion}/scatter-histogram-boxwhisker.pdf --pages . 319 -- ${caprion}/${protein}/${protein}.pdf

# net use W: \\ME-FILER1\GROUPS$\MGEU
# W:\Factors\INTERVAL\Caprion_proteomics
