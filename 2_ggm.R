#!/usr/bin/bash

options(width=200)
suppressMessages(library(Biobase))
suppressMessages(library(corpcor))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(GeneNet))
suppressMessages(library(igraph))
suppressMessages(library(RCy3))
suppressMessages(library(Rgraphviz))
suppressMessages(library(VennDiagram))

load("ZWK.rda")
load("ZYQ.rda")
load("UDP.rda")

protein_all <- Biobase::combine(protein_ZWK,protein_ZYQ) %>%
               Biobase::combine(protein_UDP)
peptide_all <- Biobase::combine(peptide_ZWK,peptide_ZYQ) %>%
               Biobase::combine(peptide_UDP)

sumstat_batch <- function(type,suffix,batch,method="ggm")
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

cor_batch <- function(n=1)
{
  z <- list()
  for (pp in c("protein","peptide")[1:n])
  {
    pcor_ZWK <- sumstat_batch(pp,"ZWK","batch1")
    pcor_ZYQ <- sumstat_batch(pp,"ZYQ","batch2")
    pcor_UDP <- sumstat_batch(pp,"UDP","batch3")
    pcor <- bind_rows(pcor_ZWK,pcor_ZYQ,pcor_UDP)
    z[[pp]] <- pcor
  }
  z
}

stats_batch <- cor_batch()
lapply(stats_batch,dim)

graph_make <- function()
{
  edges <- ggm.list.edges(pcor) %>%
           filter(node1!=node2)
  print(head(edges))
  graph <- network.make.graph(edges,labels)
  plot(graph,"fdp")
}

ggm_batch <- function(type,suffix)
{
  match.suffix <- match(suffix,c("ZWK","ZYQ","UDP"))
  pcor <- subset(stats_batch$protein,batch==switch(match.suffix,"batch1","batch2","batch3")) %>%
          select(sub("\\b(^[0-9])","\\X\\1",featureNames(get(paste(type,suffix,sep="_"))))) %>%
          as.matrix()
  labels <- colnames(pcor)
# graph_make()
  nodes <- ncol(pcor)
  pdf(file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".pdf")))
  tests <- network.test.edges(pcor)
  net <- extract.network(tests, cutoff.ggm=0.05/(nodes*(nodes-1)/2))
  graph <- network.make.graph(net,labels)
  save(graph,file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".graph")))
  g <- graph_from_graphnel(graph)
  plot(g)
  suid_corrIGraph <- createNetworkFromIgraph(g,"corrIGraph")
  layoutNetwork("attribute-circle")
  exportImage(file.path("~/Caprion/pilot/work","corrIGraph.pdf"),type="PDF",resolution=300,height=8,width=12,units="in",overwriteFile=TRUE)
  exportNetwork(file.path("~/Caprion/pilot/work","corrIGraph.sif"))
  saveSession(file.path("~/Caprion/pilot/work","corrIGraph.cys"),overwriteFile=TRUE)
  deleteNetwork(suid_corrIGraph)
  dev.off()
}

ggm_batch("protein","ZYQ")

ggm_all <- function(type,suffix)
{
  es <- get(paste(type,suffix,sep="_"))
  d_raw <- t(exprs(es))
  exclude <- names(apply(d_raw,2,sum))[is.na(apply(d_raw,2,sum))]
  d <- d_raw[, !colnames(d_raw) %in% exclude]
  p <- unclass(ggm.estimate.pcor(d))
  colnames(p) <- rownames(p) <- sub("\\b(^[0-9])","\\X\\1",colnames(d))
  labels <- colnames(p)
  nodes <- ncol(p)
  pdf(file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".pdf")))
  tests <- network.test.edges(p)
  net <- extract.network(tests, cutoff.ggm=0.05/(nodes*(nodes-1)/2))
  graph <- network.make.graph(net,labels)
  save(graph,file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".graph")))
  g <- graph_from_graphnel(graph)
  plot(g)
  dev.off()
}

ggm_all("protein","all")
