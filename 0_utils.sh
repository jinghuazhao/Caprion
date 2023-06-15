#!/usr/bin/bash

Rscript -e '
  options(width=200)
  suppressMessages(library(Biobase))
  suppressMessages(library(dplyr))
  suppressMessages(library(Hmisc))
  suppressMessages(library(pheatmap))
  filter(pQTLdata::caprion,Protein=="INHBE_HUMAN") %>%
  select(Protein,Accession,Gene,Protein.Description)
  protein_peptide <- function(protein="INHBE",suffix="ZWK")
  {
    cat(suffix,"\n")
    dir <- "~/rds/projects/Caprion_proteomics"
    load(file.path(dir,"pilot",paste(suffix,"rda",sep=".")))
    n <- paste("protein",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% paste(protein,"HUMAN",sep="_")
    prot <- matrix(p[g,],nrow=1,dimnames=list(protein,names(p[g,]))) %>%
            data.frame
    n <- paste("dr",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% paste(protein,"HUMAN",sep="_")
    dr <- matrix(p[g,],nrow=1,dimnames=list(protein,names(p[g,]))) %>%
          data.frame
    n <- paste("mapping",suffix,sep="_")
    m <- subset(get(n),grepl(protein,Protein))
    igID <- m[["Isotope.Group.ID"]]
    n <- paste("peptide",suffix,sep="_")
    p <- exprs(get(n))
    g <- rownames(p) %in% igID
    pept <- data.frame(p[g,])
    prot_pept <- bind_rows(pept,prot)
    n <- rownames(prot_pept)
    d <- t(prot_pept) %>%
         data.frame() %>%
         setNames(n)
    cat("Protein\n")
    s <- summary(lm(INHBE~d[["442593377"]]+d[["442626845"]]+d[["442628596"]]+d[["442658425"]],data=d))
    print(s,digits=3)
    dr_pept <- bind_rows(pept,dr)
    n <- rownames(dr_pept)
    d <- t(dr_pept) %>%
         data.frame() %>%
         setNames(n)
    cat("Protein_DR\n")
    s <- summary(lm(INHBE~d[["442593377"]]+d[["442626845"]]+d[["442628596"]]+d[["442658425"]],data=d))
    print(s,digits=3)
    opar <- par()
    png(file.path(dir,"analysis","work",paste(protein,suffix,"dist.png",sep="-")),
        width=12,height=10,units="in",pointsize=4,res=300)
    par(mar=c(15,10,5,5), font=2, font.lab = 5, font.axis = 5)
    source("https://raw.githubusercontent.com/jinghuazhao/tests/main/Hmisc/hist.data.frame.R")
    hist.data.frame(d,cex.axis=2.5,cex.lab=2.5,cex.mtext=2.5,cex.names=2.5,ylab=expression("Frequency"))
    dev.off()
    par(opar)
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
'

cd ~/Caprion/analysis/peptide/INHBE/
cat <(gunzip -c INHBE-?-442628596.fastGWA.gz | head -1 | paste <(echo Batch-peptide) -) \
    <(zgrep -w rs11172187 INHBE-?-442628596.fastGWA.gz) | \
sed 's/fastGWA.gz:/\t/' | \
Rscript -e '
  rs11172187 <- read.table("stdin",header=TRUE)
  knitr::kable(rs11172187,caption="Effect sizes of rs11172187",digits=3)
'
cd -
