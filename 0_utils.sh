#!/usr/bin/bash

function pgwas()
{
  export protein=${1}
  export pqtl=${2}
  cd ~/Caprion/analysis/pgwas
  cat <(gunzip -c caprion-?-${protein}.fastGWA.gz | head -1 | paste <(echo Batch-Protein) -|sed 's/-/\t/') \
      <(zgrep -w ${pqtl} caprion-?-${protein}.fastGWA.gz) | \
  sed 's/caprion-//;s/.fastGWA.gz:/\t/;s/-/\t/' | \
  Rscript -e '
    suppressMessages(library(dplyr))
    pqtl <- Sys.getenv("pqtl")
    d <- read.table("stdin",header=TRUE) %>%
         arrange(Batch)
    knitr::kable(d,caption=paste("Effect sizes of",pqtl),digits=3)
  '
  cd ~/Caprion/analysis/peptide/${protein}/
  cat <(gunzip -c INHBE-?-*.fastGWA.gz | head -1 | paste <(echo Batch-Peptide) -|sed 's/-/\t/') \
      <(zgrep -w ${pqtl} *fastGWA.gz) | \
  sed 's/INHBE-//;s/.fastGWA.gz:/\t/;s/-/\t/' | \
  Rscript -e '
    suppressMessages(library(dplyr))
    pqtl <- Sys.getenv("pqtl")
    d <- read.table("stdin",header=TRUE) %>%
         arrange(Peptide,Batch)
    knitr::kable(d,caption=paste("Effect sizes of",pqtl),digits=3)
  '
}

function hist_corr_lm()
{
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
    cat("\n**",suffix,"**\n",sep="")
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
    cat("\nProtein")
    s <- summary(lm(INHBE~d[["442593377"]]+d[["442626845"]]+d[["442628596"]]+d[["442658425"]],data=d))
    print(knitr::kable(coef(s),digits=3))
    dr_pept <- bind_rows(pept,dr)
    n <- rownames(dr_pept)
    d <- t(dr_pept) %>%
         data.frame() %>%
         setNames(n)
    cat("\nProtein_DR")
    s <- summary(lm(INHBE~d[["442593377"]]+d[["442626845"]]+d[["442628596"]]+d[["442658425"]],data=d))
    print(knitr::kable(coef(s),digits=3))
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
}

pgwas INHBE rs149830883
pgwas INHBE rs11172187
