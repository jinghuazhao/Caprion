module load gcc/5
R --no-save <<END
  rmarkdown::render("ae_test.Rmd", c("pdf_document","html_document"))
END
