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
suppressMessages(library(quantro))
suppressMessages(library(sva))
suppressMessages(library(VennDiagram))
suppressMessages(library(visNetwork))

load("~/Caprion/pilot/ZWK.rda")
load("~/Caprion/pilot/ZYQ.rda")
load("~/Caprion/pilot/UDP.rda")
load("~/Caprion/pilot/UHZ.rda")

protein_all <- Biobase::combine(protein_ZWK,protein_ZYQ) %>%
               Biobase::combine(protein_UDP)
peptide_all <- Biobase::combine(peptide_ZWK,peptide_ZYQ) %>%
               Biobase::combine(peptide_UDP)
protein_dr_all <- Biobase::combine(dr_ZWK,dr_ZYQ) %>%
                  Biobase::combine(dr_UDP)
save(protein_all,protein_dr_all,peptide_all,file="~/Caprion/pilot/es.rda")

protein_all <- Biobase::combine(protein_ZWK,protein_ZYQ) %>%
               Biobase::combine(protein_UDP) %>%
               Biobase::combine(protein_UHZ)
peptide_all <- Biobase::combine(peptide_ZWK,peptide_ZYQ) %>%
               Biobase::combine(peptide_UDP) %>%
               Biobase::combine(peptide_UHZ)
protein_dr_all <- Biobase::combine(dr_ZWK,dr_ZYQ) %>%
                  Biobase::combine(dr_UDP) %>%
                  Biobase::combine(dr_UHZ)
save(protein_all,protein_dr_all,peptide_all,file="~/Caprion/analysis/work/es.rda")

load("~/Caprion/pilot/work/es.rda")
prot3 <- subset(protein_all,!featureNames(protein_all)%in%featureNames(protein_UHZ))
ZYQ.na <- paste(c("BROX","CT027","GHRL","PSB6"),"_HUMAN",sep="")
UDP.na <- paste(c("BROX","NCF2","SEM7A"),"_HUMAN",sep="")
comm_all <- subset(protein_all,!featureNames(protein_all)%in%featureNames(prot3))
load("~/Caprion/analysis/work/es.rda")
prot4 <- subset(protein_all,featureNames(protein_all)%in%featureNames(comm_all) &
                           !featureNames(protein_all)%in%c(ZYQ.na,UDP.na))
edata <- prot4[,!is.na(apply(exprs(prot4),2,sum))]
coln <- colnames(edata)
coln_ZWK <- grepl("ZWK",coln)
coln_ZYQ <- grepl("ZYQ",coln)
coln_UDP <- grepl("UDP",coln)
coln_UHZ <- grepl("UHZ",coln)
coln[coln_ZWK] <- 1
coln[coln_ZYQ] <- 2
coln[coln_UDP] <- 3
coln[coln_UHZ] <- 4

png("~/work/matboxplot.png",width=10,height=8,units="in",res=300)
par(cex=0.8,pch=19)
matboxplot(edata,groupFactor=coln, ylab="Protein measurement")
dev.off()
png("~/work/matdensity.png",width=10,height=8,units="in",res=300)
matdensity(edata, coln, xlab = " ", ylab = "density", ylim=c(0,2),
           main = "Protein levels", brewer.n = 8, brewer.name = "Dark2")
legend('top', c("ZWK", "ZYQ", "UDP", "UHZ"), col = c(1, 2, 3, 4), lty = 1, lwd = 3)
dev.off()
png("~/work/combat.png",width=10,height=8,units="in",res=300)
combat_edata <- ComBat(dat=edata, batch=coln, mod=NULL, par.prior=TRUE, prior.plots=TRUE, ref.batch=3)
dev.off()

quantro_qtest  <- function(dat)
{
   qtest <- quantro(dat,groupFactor=coln)
   print(anova(qtest))
   print(quantroStat(qtest))
   qsp <- quantroStatPerm(qtest)
   print(quantroPvalPerm(qtest))
   matboxplot(dat, groupFactor=coln, cex=0.4, pch=19, xaxt="n",
              main="Boxplots of all proteins", xlab="Sample", ylab="Abundance level",
              brewer.n=8, brewer.name="Dark2")
   legend('ZWK', c("ZYQ", "UDP","UHZ"), col = 1:4, lty = 1, lwd = 3)
   q2q3 <- apply(dat,2,quantile,c(0.25,0.75))
   outside <- dat<q2q3[1,]|dat>q2q3[2,]
   outside_prot <- apply(outside+0,1,sum,na.rm=TRUE)
   outside_id <- apply(outside+0,2,sum,na.rm=TRUE)
   upper <- dat>q2q3[2,]
   upper_prot <- apply(upper+0,1,sum,na.rm=TRUE)
   upper_id <- apply(upper+0,2,sum,na.rm=TRUE)
   lower <- dat<q2q3[1,]
   lower_prot <- apply(lower+0,1,sum,na.rm=TRUE)
   lower_id <- apply(lower+0,2,sum,na.rm=TRUE)
   x <- list(outside,outside_prot,outside_id,upper,upper_prot,upper_id,lower,lower_prot,lower_id)
   setNames(x,c("outside","outside_prot","outside_id","upper","upper_prot","upper_id","lower","lower_prot","lower_id"))
}

png("~/work/quantro.png",width=10,height=8,units="in",res=300)
par(mfrow=c(2,1))
z <- quantro_qtest(edata)
n <- quantro_qtest(combat_edata)
dev.off()

fcheck <- function(es)
{
  fn <- featureNames(es)
  d <- exprs(es)
  rs <- apply(d,1,sum)
  nna <- names(rs[is.na(rs)])
  list(fn=fn,nna=nna,band=d[nna,])
}

all <- fcheck(protein_all)
all$fn[grepl("PP2|^6|^1",all$fn)]
zwk <- fcheck(protein_ZWK)
zyq <- fcheck(protein_ZYQ)
udp <- fcheck(protein_UDP)
lapply(list(all$fn,zwk$fn,zyq$fn,udp$fn),length)
setdiff(all$fn,zwk$fn)
setdiff(zwk$fn,zyq$fn)
setdiff(zwk$fn,udp$fn)
setdiff(zyq$fn,udp$fn)

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
  q4 <- with(net,quantile(pcor))
  c4 <- with(net,cut(pcor,q4))
  edges <- with(net,data.frame(from=node1,to=node2,value=30*pcor,
                               color=c("#0000FF","#9999FF","#00FF00","#FF9999","#FF0000")[c4]))
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
  visSave(network,file=file.path("~/Caprion/analysis/work","protein_all.html"),selfcontained=TRUE)
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
  pdf(file=file.path("~/Caprion/analysis/work",paste0(type,"_",suffix,".pdf")))
  tests <- network.test.edges(pcor)
  net <- extract.network(tests, cutoff.ggm=0.05/(nodes*(nodes-1)/2))
  graph <- network.make.graph(net,labels)
  g <- graph_from_graphnel(graph)
  save(tests,net,graph,g,file=file.path("~/Caprion/analysis/work",paste0(type,"_",suffix,".graph")))
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
  exportImage(file.path("~/Caprion/analysis/work","corrIGraph.pdf"),type="PDF",resolution=300,height=8,width=12,units="in",overwriteFile=TRUE)
  exportNetwork(file.path("~/Caprion/analysis/work","corrIGraph.sif"))
  saveSession(file.path("~/Caprion/analysis/work","corrIGraph.cys"),overwriteFile=TRUE)
  deleteNetwork(suid_corrIGraph)
}
