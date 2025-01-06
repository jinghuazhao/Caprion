#!/usr/bin/bash

# https://platform.opentargets.org/api
# ENSG00000164308

module load ceuadmin/R
Rscript -e '
  ERAP2 <- subset(pQTLdata::caprion,Gene=="ERAP2")
  ERAP2$ensGenes
  data <- jsonlite::fromJSON("ERAP2.json")
  diseases_data <- data$data$target$associatedDiseases$rows
  diseases_data <- tidyr::unnest(diseases_data, cols = everything(), names_sep = "_")
  write.table(diseases_data, file = "ERAP2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
'

# https://platform-docs.opentargets.org/data-access/graphql-api

## Python: Used directly without change

python3 <<END
import requests
import json

gene_id = "ENSG00000164308"
query_string = """
  query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      biotype
      geneticConstraint {
        constraintType
        exp
        obs
        score
        oe
        oeLower
        oeUpper
      }
      tractability {
        label
        modality
        value
      }
    }
  }
"""
variables = {"ensemblId": gene_id}
base_url = "https://api.platform.opentargets.org/api/v4/graphql"
r = requests.post(base_url, json={"query": query_string, "variables": variables})
print(r.status_code)
api_response = json.loads(r.text)
print(api_response)
END

## R: necessary to get around httr::content(r) for iconvlist() with iconv().

Rscript -e '
library(httr)
library(jsonlite)

gene_id <- "ENSG00000164308"
query_string = "
  query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      biotype
      geneticConstraint {
        constraintType
        exp
        obs
        score
        oe
        oeLower
        oeUpper
      }
      tractability {
        label
        modality
        value
      }
    }
  }
"
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
variables <- list("ensemblId" = gene_id)
post_body <- list(query = query_string, variables = variables)
r <- httr::POST(url=base_url, body=post_body, encode='json')
data <- iconv(r, "latin1", "ASCII")
content <- jsonlite::fromJSON(data)
'
