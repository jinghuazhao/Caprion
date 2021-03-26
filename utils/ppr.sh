#!/usr/bin/bash

Rscript -e "knitr::knit('ppr.Rmd')"
pandoc ppr.md -o ppr.docx
