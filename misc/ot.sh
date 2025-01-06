#!/usr/bin/bash

# https://platform.opentargets.org/api
# ENSG00000164308

module load ceuadmin/R
Rscript -e '
  ERAP2 <- subset(pQTLdata::caprion,Gene=="ERAP2")
  ERAP2$ensGenes
  data <- jsonlite::fromJSON("ERAP2.json")
  diseases_data <- data$data$target$associatedDiseases$rows
  library(tidyr)
  diseases_data <- unnest(diseases_data, cols = everything(), names_sep = "_")
  write.table(diseases_data, file = "ERAP2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
'
