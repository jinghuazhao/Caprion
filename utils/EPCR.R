library(pQTLtools)

annotation <- function()
## Annotations
{
  pag <- subset(caprion,Accession=="Q9UNN8")
#    Protein Accession  Gene                                                          Protein.Description
# EPCR_HUMAN    Q9UNN8 PROCR Endothelial protein C receptor (Activated protein C receptor) (APC receptor)

  subset(hg19Tables,hgncSym=="PROCR")
# chrom chromStart chromEnd strand    acc uniprotName geneName geneSynonyms hgncSym         ensGene
# chr20   33759957 33764613      + Q9UNN8  EPCR_HUMAN    PROCR         EPCR   PROCR ENSG00000101000

  qqman <- pag[c("Accession","Gene")]
  write.table(qqman,file="qqman.list",quote=FALSE,row.names=FALSE,col.names=FALSE)
}

load("caprion.rda")

extract <- function(prots=c("EPCR_HUMAN","PROC_HUMAN"))
{
# Isotope.Group.ID, Modified.Peptide.Sequence
  d1 <- subset(Normalized_Peptides,Protein %in% prots)
  piptides <- t(d1[,-(1:6)])
  colnames(piptides) <- paste0(gsub("HUMAN","",d1[,2]),d1[,1])
# Protein
  d2 <- subset(Protein_All_Peptides,Protein %in% prots)
  all <- t(d2[,-1])
  colnames(all) <- gsub("HUMAN","All",d2[,1])
# Protein_DR_Filt
  d3 <- subset(Protein_All_Peptides,Protein %in% prots)
  dr <- t(d3[,-1])
  colnames(dr) <- gsub("HUMAN","DR",d3[,1])
  cbind(piptides,all,dr)
}

epcr_proc <- extract()

