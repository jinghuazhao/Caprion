library(pQTLtools)

pag <- subset(caprion,Accession=="Q9UNN8")
#    Protein Accession  Gene                                                          Protein.Description
# EPCR_HUMAN    Q9UNN8 PROCR Endothelial protein C receptor (Activated protein C receptor) (APC receptor)

subset(hg19Tables,hgncSym=="PROCR")
# chrom chromStart chromEnd strand    acc uniprotName geneName geneSynonyms hgncSym         ensGene
# chr20   33759957 33764613      + Q9UNN8  EPCR_HUMAN    PROCR         EPCR   PROCR ENSG00000101000

qqman <- pag[c("Accession","Gene")]
write.table(qqman,file="qqman.list",quote=FALSE,row.names=FALSE,col.names=FALSE)
