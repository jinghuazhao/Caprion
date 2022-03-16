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
suppressMessages(library(visNetwork))

load("~/Caprion/pilot/ZWK.rda")
load("~/Caprion/pilot/ZYQ.rda")
load("~/Caprion/pilot/UDP.rda")

protein_all <- Biobase::combine(protein_ZWK,protein_ZYQ) %>%
               Biobase::combine(protein_UDP)
peptide_all <- Biobase::combine(peptide_ZWK,peptide_ZYQ) %>%
               Biobase::combine(peptide_UDP)

ggm_all <- function(type,suffix)
{
  es <- get(paste(type,suffix,sep="_"))
  d_raw <- t(exprs(es))
  d_sum <- apply(d_raw,2,sum)
  exclude <- names(d_sum)[is.na(d_sum)]
  d <- d_raw[, !colnames(d_raw) %in% exclude]
  p <- unclass(ggm.estimate.pcor(d))
  colnames(p) <- rownames(p) <- gsub("_HUMAN","",colnames(d))
  labels <- colnames(p)
  nnodes <- ncol(p)
  pdf(file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".pdf")))
  tests <- network.test.edges(p)
  e <- extract.network(tests, cutoff.ggm=0.05/(nnodes*(nnodes-1)/2))
  id=sort(unique(c(e[["node1"]],e[["node2"]])))
  net <- mutate(e,label1=labels[node1],label2=labels[node2])
  graph <- network.make.graph(net,labels)
  g <- graph_from_graphnel(graph)
  save(tests,net,graph,g,file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".graph")))
  plot(g)
  dev.off()
  title <- list(text="Gaussian graphical models of proteins",
                style="font-family:Arial;color:black;font-size:30px;text-align:center;")
  nodes <- data.frame(id,label=labels[id],shape="box")
  edges <- with(net,data.frame(from=node1,to=node2,value=30*pcor))
  nodesId <- list(enabled = TRUE,
                  style='width: 200px; height: 26px;
                         background: #f8f8f8;
                         color: darkblue;
                         border:none;
                         outline:none;')
  network <- visNetwork(nodes,edges,width=1500,height=1250,main=title) %>%
             visOptions(highlightNearest=TRUE, nodesIdSelection=nodesId) %>%
             visInteraction(navigationButtons=TRUE) %>%
             visIgraphLayout(type="full") %>%
             visNodes(size=30)
  visSave(network,file=file.path("~/Caprion/pilot/work","protein_all.html"),selfcontained=TRUE)
}

ggm_all("protein","all")

sumstat_batch <- function(type,suffix,method="ggm")
{
  es <- get(paste(type,suffix,sep="_"))
  cat(type,suffix,"No. features =",length(featureNames(es)),"\n")
  d <- t(exprs(es))
  match.type <- match(method,c("cor","cor.shrink","ggm"))
  p <- switch(match.type,
              {estimate.lambda(d);cor2pcor(cor(d,use="everything",method="pearson"))},
              cor2pcor(cor.shrink(d)),
              unclass(ggm.estimate.pcor(d)))
  colnames(p) <- rownames(p) <- featureNames(es)
  p
}

ggm_batch <- function(type,suffix)
{
  pcor <- sumstat_batch(type,suffix)
  labels <- colnames(pcor)
# graph_make()
  nodes <- ncol(pcor)
  pdf(file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".pdf")))
  tests <- network.test.edges(pcor)
  net <- extract.network(tests, cutoff.ggm=0.05/(nodes*(nodes-1)/2))
  graph <- network.make.graph(net,labels)
  g <- graph_from_graphnel(graph)
  save(tests,net,graph,g,file=file.path("~/Caprion/pilot/work",paste0(type,"_",suffix,".graph")))
  plot(g)
  dev.off()
# graph_RCy3()
}

ggm_batch("protein","ZYQ")
graph_make <- function()
{
  edges <- ggm.list.edges(pcor) %>%
           filter(node1!=node2)
  print(head(edges))
  graph <- network.make.graph(edges,labels)
  plot(graph,"fdp")
}

graph_RCy3 <- function()
{
  suid_corrIGraph <- createNetworkFromIgraph(g,"corrIGraph")
  layoutNetwork("attribute-circle")
  exportImage(file.path("~/Caprion/pilot/work","corrIGraph.pdf"),type="PDF",resolution=300,height=8,width=12,units="in",overwriteFile=TRUE)
  exportNetwork(file.path("~/Caprion/pilot/work","corrIGraph.sif"))
  saveSession(file.path("~/Caprion/pilot/work","corrIGraph.cys"),overwriteFile=TRUE)
  deleteNetwork(suid_corrIGraph)
}
