#!/usr/bin/bash

Rscript -e "knitr::knit('ppr.Rmd')"
pandoc ppr.md -o ppr.docx

# Rscript -e 'rmarkdown::render("ppr.Rmd", output_format = "html_document")'
# Rscript -e 'rmarkdown::render("ppr.Rmd", output_format = "pdf_document")'
# Rscript -e 'rmarkdown::render("ppr.Rmd", output_format = "word_document")'
