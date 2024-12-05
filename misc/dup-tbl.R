#!/usr/local/Cluster-Apps/ceuadmin/R/4.4.2-icelake/bin/Rscript

suppressMessages(library(dplyr))
analysis <- Sys.getenv("analysis")
lst <- read.delim(file.path(analysis,"dup","dup.list"),
                  col.names=c("Isotope.Group.ID","Protein","Modified.Peptide.Sequence","Monoisotopic","Time","Charge"))
m <- read.delim(file.path(analysis,"dup","dup.merge"))
s <- read.delim(file.path(analysis,"dup","dup.signals")) %>%
     dplyr::select(trait,P,SNP)
signals <- left_join(m,s)
tbl <- dplyr::left_join(lst,signals,by=c("Isotope.Group.ID"="trait"))
write.table(tbl, file = file.path(analysis,"dup","dup.tbl"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
