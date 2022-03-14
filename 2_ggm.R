#!/usr/bin/bash

options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(arrayQualityMetrics))
suppressMessages(library(corpcor))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(VennDiagram))

load("ZWK.rda")
load("ZYQ.rda")
load("UDP.rda")

sumstats <- function(type,suffix,batch,method="cor.shrink")
{
  es <- get(paste(type,suffix,sep="_"))
  cat(type,batch,"No. features =",length(featureNames(es)),"\n")
  d <- t(exprs(es))
  match.type <- match(method,c("cor","cor.shrink"))
  r <- switch(match.type,{estimate.lambda(d);cor(d,use="everything",method="pearson")},cor.shrink(d))
  p <- cor2pcor(r)
  colnames(p) <- colnames(r)
  rownames(p) <- rownames(r)
  data.frame(p,batch=batch)
}

cors <- function(n=1)
{
  z <- list()
  for (pp in c("protein","peptide")[1:n])
  {
    pcor_ZWK <- sumstats(pp,"ZWK","batch1")
    pcor_ZYQ <- sumstats(pp,"ZYQ","batch2")
    pcor_UDP <- sumstats(pp,"UDP","batch3")
    pcor <- bind_rows(pcor_ZWK,pcor_ZYQ,pcor_UDP)
    z[[pp]] <- pcor
  }
  z
}
stats <- cors()
lapply(stats,dim)

suppressMessages(library("GeneNet"))
suppressMessages(library("Rgraphviz"))
pcor_ZWK <- subset(stats$protein,batch=="batch1") %>%
            select(sub("\\b(^[0-9])","\\X\\1",featureNames(protein_ZWK))) %>%
            as.matrix()
label_ZWK <- colnames(pcor_ZWK)
edges_ZWK <- ggm.list.edges(pcor_ZWK)
graph_ZWK <- network.make.graph(edges_ZWK,label_ZWK)
plot(graph_ZWK,"fdp")
