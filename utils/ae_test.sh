module load gcc/5
R --no-save <<END
rmarkdown::render("ae_test.Rmd", c("html_document", "pdf_document"))
END

