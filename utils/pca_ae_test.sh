R --no-save <<END
rmarkdown::render("pca_ae_test.Rmd", c("html_document", "pdf_document"))
END

