#!/usr/bin/bash

Rscript -e '
  suppressMessages(library(Biobase))
  suppressMessages(library(dplyr))
  suppressMessages(library(Hmisc))
  suppressMessages(library(pheatmap))
  protein_peptide <- function(protein="INHBE",suffix="ZWK")
  {
    dir <- "~/rds/projects/Caprion_proteomics"
    load(file.path(dir,"pilot",paste(suffix,"rda",sep=".")))
    n <- paste("protein",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% paste(protein,"HUMAN",sep="_")
    prot <- matrix(p[g,],nrow=1,dimnames=list(protein,names(p[g,]))) %>%
            data.frame
    n <- paste("mapping",suffix,sep="_")
    m <- subset(get(n),grepl(protein,Protein))
    igID <- m[["Isotope.Group.ID"]]
    n <- paste("peptide",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% igID
    pept <- data.frame(p[g,])
    prot_pept <- bind_rows(pept,prot)
    png(file.path(dir,"analysis","work",paste(protein,suffix,"dist.png",sep="-")),
        width=12,height=10,units="in",pointsize=4,res=300)
    par(mar=c(10,5,3,3), font.lab = 3.5, font.axis = 3.5)
    n <- rownames(prot_pept)
    d <- t(prot_pept) %>%
         data.frame() %>%
         setNames(n)
    hist.data.frame(d,cex.axis=2.5,cex.lab=2.5,cex.mtext=2.5,cex.names=2.5)
    dev.off()
    png(file.path(dir,"analysis","work",paste(protein,suffix,"corr.png",sep="-")),
        width=12,height=10,units="in",pointsize=4,res=300)
    pheatmap(cor(t(prot_pept)),display_numbers=TRUE,fontsize=24)
    dev.off()
    write.csv(m[c("Isotope.Group.ID","Modified.Peptide.Sequence","Protein")],
              file=file.path(dir,"analysis","work",paste(protein,suffix,"mapping.csv",sep="-")),
              quote=FALSE,row.names=FALSE)
    prot_pept
  }
  zwk <- protein_peptide()
  zyq <- protein_peptide(suffix="ZYQ")
  udp <- protein_peptide(suffix="UDP")
  filter(pQTLdata::caprion,Protein=="INHBE_HUMAN") %>%
  select(Protein,Accession,Gene,Protein.Description)
'
