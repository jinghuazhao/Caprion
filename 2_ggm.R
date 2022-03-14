#!/usr/bin/bash

options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(arrayQualityMetrics))
suppressMessages(library(corpcor))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(GeneNet))
suppressMessages(library(Rgraphviz))
suppressMessages(library(VennDiagram))

load("ZWK.rda")
load("ZYQ.rda")
load("UDP.rda")

sumstats <- function(type,suffix,batch,method="ggm")
{
  es <- get(paste(type,suffix,sep="_"))
  cat(type,batch,"No. features =",length(featureNames(es)),"\n")
  d <- t(exprs(es))
  match.type <- match(method,c("cor","cor.shrink","ggm"))
  p <- switch(match.type,
              {estimate.lambda(d);cor2pcor(cor(d,use="everything",method="pearson"))},
              cor2pcor(cor.shrink(d)),
              unclass(ggm.estimate.pcor(d)))
  colnames(p) <- rownames(p) <- sub("\\b(^[0-9])","\\X\\1",featureNames(es))
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

graph_all <- function()
{
  edges <- ggm.list.edges(pcor) %>%
           filter(node1!=node2)
  print(head(edges))
  graph <- network.make.graph(edges,labels)
  plot(graph,"fdp")
}

ggm <- function(type,suffix)
{
  match.suffix <- match(suffix,c("ZWK","ZYQ","UDP"))
  pcor <- subset(stats$protein,batch==switch(match.suffix,"batch1","batch2","batch3")) %>%
          select(sub("\\b(^[0-9])","\\X\\1",featureNames(get(paste(type,suffix,sep="_"))))) %>%
          as.matrix()
  labels <- colnames(pcor)
# graph_all()
  nodes <- ncol(pcor)
  tests <- network.test.edges(pcor)
  net <- extract.network(tests, cutoff.ggm=0.05/(nodes*(nodes-1)/2))
  graph <- network.make.graph(net,labels)
  plot(graph,"fdp")
}

ggm("protein","ZWK")
