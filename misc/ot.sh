#!/usr/bin/bash

# https://platform.opentargets.org/api
# replaced with "ENSG00000164308"

module load ceuadmin/R
# Saved output from the web
Rscript -e '
  ERAP2 <- subset(pQTLdata::caprion,grepl("ERAP2",Gene))
  ERAP2$ensGenes
  data <- jsonlite::fromJSON("ERAP2.json")
  diseases_data <- data$data$target$associatedDiseases$rows
  diseases_data <- tidyr::unnest(diseases_data, cols = everything(), names_sep = "_")
  write.table(diseases_data, file = "ERAP2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
'
# Replicate from scratch
Rscript -e '
library(httr)
library(jsonlite)
gene_id <- "ENSG00000164308"
query_string <- "
  query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      associatedDiseases {
        count
        rows {
          disease {
            id
            name
          }
          datasourceScores {
            id
            score
          }
        }
      }
    }
  }
"
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
variables <- list("ensemblId" = gene_id)
post_body <- list(query = query_string, variables = variables)
r <- httr::POST(url = base_url, body = post_body, encode = "json")
if (status_code(r) == 200) {
  data <- iconv(r, "latin1", "ASCII")
  content <- jsonlite::fromJSON(data)
  write(data,file="ERAP2.json")
} else {
  print(paste("Request failed with status code", status_code(r)))
}
# multiple genes
library(httr)
library(jsonlite)
benchmarks <- subset(pQTLdata::caprion,grepl("A1BG|APOB|ERAP2|PROC",Gene))
subset(benchmarks,select=c(Gene,ensGenes))
gene_ids <- benchmarks$ensGenes

query_string <- "
  query target($ensemblIds: [String!]!){
    targets(ensemblIds: $ensemblIds){
      id
      approvedSymbol
      associatedDiseases {
        count
        rows {
          disease {
            id
            name
          }
          datasourceScores {
            id
            score
          }
        }
      }
    }
  }
"
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
variables <- list("ensemblIds" = gene_ids)
post_body <- list(query = query_string, variables = variables)
r <- httr::POST(url = base_url, body = post_body, encode = "json")
data <- iconv(r, "latin1", "ASCII")
content <- jsonlite::fromJSON(data)
'

# https://platform-docs.opentargets.org/data-access/graphql-api

## Python: Used directly without change
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
r <- httr::POST(url=base_url, body=post_body, encode="json")
data <- iconv(r, "latin1", "ASCII")
content <- jsonlite::fromJSON(data)
'
