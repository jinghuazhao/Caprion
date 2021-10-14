UDP <- function()
{
  r <- sapply(1:length(featureNames(protein_UDP)),function(r) {
              protein_UDP_r <- protein_UDP[r,]
              fn <- paste0("invnormal(",sub("(^[0-9])","X\\1",featureNames(protein_UDP_r)),")")
              f <- paste(fn,"~ agePulse + sexPulse + classification")
              z <- lm(as.formula(f),data=protein_UDP_r, na.action=na.exclude)
              resid(z)
            })
  colnames(r) <- featureNames(protein_UDP)
  d <- data.frame(r)
  d <- d %>%
       mutate(caprion_id=rownames(r)) %>%
       left_join(pData(protein_UDP)[c("caprion_id","Affymetrix_gwasQC_bl")]) %>%
       select(Affymetrix_gwasQC_bl,caprion_id,setdiff(names(d),c("Affymetrix_gwasQC_bl","caprion_id"))) %>%
       filter(!caprion_id %in%c("UDP0138","UDP0481"))
  names(d) <- c("FID","IID",featureNames(protein_UDP))
  write.table(d,file=file.path(caprion,"data3","UDP.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
}

caprion <- Sys.getenv("caprion")
library(gap)
load("UDP.rda")
