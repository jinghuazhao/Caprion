!/usr/bin/bash

export analysis=~/Caprion/analysis

function tables()
{
R --no-save <<END
require(openxlsx)
suppressMessages(library(dplyr))
suppressMessages(library(gap))
require(stringr)
options(width=2000)
options("openxlsx.borderColour"="#4F80BD")
analysis <- Sys.getenv("analysis")
freq <- read.table("~/Caprion/analysis/bgen/caprion.freq",col.names=c("SNP","REF","ALT","ALT_FREQS"))
tbl <- read.delim(file.path(analysis,"work","tbl.tsv")) %>%
       dplyr::select(prot,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,Direction,HetISq,logHetP,N) %>%
       dplyr::mutate(Allele1=toupper(Allele1),Allele2=toupper(Allele2)) %>%
       dplyr::rename(Protein=prot,SNP=MarkerName,EA=Allele1,OA=Allele2,EAF=Freq1,Log10P=log.P.)
tbl_dr <- read.delim(file.path(analysis,"work","tbl_dr.tsv")) %>%
          dplyr::select(prot,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,Direction,HetISq,logHetP,N) %>%
          dplyr::mutate(Allele1=toupper(Allele1),Allele2=toupper(Allele2)) %>%
          dplyr::rename(Protein=prot,SNP=MarkerName,EA=Allele1,OA=Allele2,EAF=Freq1,Log10P=log.P.)
Protein <- read.csv(file.path(analysis,"work","caprion.cis.vs.trans")) %>%
           dplyr::mutate(geneRegion=paste(geneChrom,":",geneStart,"-",geneEnd)) %>%
           dplyr::left_join(freq) %>%
           dplyr::arrange(prot,SNPChrom,SNPPos) %>%
           dplyr::rename(Protein=prot) %>%
           dplyr::select(-log10p,-geneChrom,-geneStart,-geneEnd,-SNPChrom,-SNPPos,-cis) %>%
           dplyr::left_join(tbl) %>%
           dplyr::select(Gene,Protein,geneRegion,SNP,Type,EA,OA,EAF,Effect,StdErr,Log10P,Direction,HetISq,logHetP,N,REF,ALT,ALT_FREQS)
Protein_dr <- read.csv(file.path(analysis,"work","caprion_dr.cis.vs.trans")) %>%
           dplyr::mutate(geneRegion=paste(geneChrom,":",geneStart,"-",geneEnd)) %>%
           dplyr::left_join(freq) %>%
           dplyr::arrange(prot,SNPChrom,SNPPos) %>%
           dplyr::rename(Protein=prot) %>%
           dplyr::select(-log10p,-geneChrom,-geneStart,-geneEnd,-SNPChrom,-SNPPos,-cis) %>%
           dplyr::left_join(tbl_dr) %>%
           dplyr::select(Gene,Protein,geneRegion,SNP,Type,EA,OA,EAF,Effect,StdErr,Log10P,Direction,HetISq,logHetP,N,REF,ALT,ALT_FREQS)
proteins <- dir("~/Caprion/analysis/peptide") # c("A1BG","APOB","EPCR","ERAP2","PROC")
Peptides <- data.frame()
for (peptide in proteins)
{
   f <- file.path(analysis,"peptide",peptide,paste0(peptide,".cis.vs.trans"))
   if(file.exists(f))
   {
     p <- read.csv(f)
     if(nrow(p)>0) Peptides <- rbind(Peptides,p)
   }
}
repl <- with(Peptides,grepl("1433$",prot))
if (nrow(Peptides[repl,])>0) Peptides[repl,"prot"] <- "1433E"
any(duplicated(freq$SNP))
dup_freq <- freq[duplicated(freq$SNP) | duplicated(freq$SNP, fromLast = TRUE), ]
dim(dup_freq)
freq_unique <- freq %>%
               dplyr::distinct(SNP, .keep_all = TRUE)
Peptides <- dplyr::left_join(Peptides,freq_unique,by="SNP") %>%
            dplyr::arrange(prot,isotope,SNPChrom,SNPPos) %>%
            dplyr::rename(Protein=prot)
snplist_0.01 <- dplyr::filter(freq,ALT_FREQS>=0.01 & ALT_FREQS<=0.99) %>%
                dplyr::pull(SNP)
Protein_0.01 <- dplyr::filter(Protein,SNP %in% snplist_0.01)
Protein_dr_0.01 <- dplyr::filter(Protein_dr,SNP %in% snplist_0.01)
Peptides_0.01 <- dplyr::filter(Peptides,SNP %in% snplist_0.01)
snplist_0.01_0.001 <- dplyr::filter(freq,(ALT_FREQS>=0.001 & ALT_FREQS<0.01)|(ALT_FREQS>0.99 & ALT_FREQS<=0.999)) %>%
                      dplyr::pull(SNP)
Protein_0.01_0.001 <- dplyr::filter(Protein,SNP %in% snplist_0.01_0.001)
Protein_dr_0.01_0.001 <- dplyr::filter(Protein_dr,SNP %in% snplist_0.01_0.001)
Peptides_0.01_0.001 <- dplyr::filter(Peptides,SNP %in% snplist_0.01_0.001)
caprion <- dplyr::mutate(Protein, ProteinSNP=paste0(Protein,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
           dplyr::rename(prot=Protein,rsid=SNP,uniprot=Gene)
caprion_dr <- dplyr::mutate(Protein_dr, ProteinSNP=paste0(Protein,"-",SNP),inv_chr_pos_a1_a2(SNP),pos=as.integer(pos)) %>%
              dplyr::rename(prot=Protein,rsid=SNP,uniprot=Gene)
same <- intersect(caprion[["ProteinSNP"]],caprion_dr[["ProteinSNP"]])
cat("Total # of chromosome X signals =",length(same[grepl("X:",same)]),"\n")
r <- cor(dplyr::filter(caprion,ProteinSNP %in% same)["Effect"],
         dplyr::filter(caprion_dr,ProteinSNP %in% same)["Effect"])
png(file.path(analysis,"work","un_dr_overlap.png"),height=8,width=8,units="in",res=300)
plot(dplyr::filter(caprion,ProteinSNP %in% same)[["Effect"]],
     dplyr::filter(caprion_dr,ProteinSNP %in% same)[["Effect"]],pch=19,
     main=paste0("Effect size comparison (N. pQTLs=",length(same),",r=",trunc(r*1e5,5)/1e5,")"),xlab="Beta",ylab="Beta_dr")
dev.off()
INF <- "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF"
plink <- "/rds/user/jhz22/hpc-work/bin/plink"
b <- list()
for(i in unique(pull(caprion,chr)))
{
   k <- dplyr::filter(caprion,chr %in% i) %>%
        dplyr::select(chr,pos,uniprot,rsid,prot)
   dr <- dplyr::filter(caprion_dr,chr %in% i) %>%
         dplyr::select(chr,pos,uniprot,rsid,prot)
   bfile <- file.path(INF,"INTERVAL","per_chr",paste0("snpid",i))
   b[[i]] <- pQTLtools::novelty_check(k,dr,ldops=list(bfile=bfile,plink=plink))
}
b[["23"]] <- mutate(b[["X"]],known.seqnames="23",query.seqnames="23")
replication <- dplyr::filter(bind_rows(b[-which(names(b)=="X")]),r2>=0.8) %>%
               dplyr::rename(known.gene=known.uniprot,query.gene=query.uniprot)
Replication <- dplyr::mutate(replication,seqnames=as.integer(known.seqnames),pos=as.integer(known.pos)) %>%
               dplyr::arrange(seqnames,pos) %>%
               dplyr::select(-known.start,-known.end,-query.seqnames,-query.start,-query.end,-seqnames,-pos)
deCODE <- read.delim(file.path(analysis,"deCODE","deCODE.tsv"))
UKB_PPP <- read.delim(file.path(analysis,"UKB_PPP","UKB_PPP.tsv"))
hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
xlsx <- "https://jhz22.user.srcf.net/Caprion/results.xlsx"
xlsx <- file.path(analysis,"reports","Supplementary-Tables.xlsx")
wb <- createWorkbook(xlsx)
addWorksheet(wb,"Summary",zoom=150)
writeData(wb,"Summary","Summary",xy=c(1,1),headerStyle=createStyle(textDecoration="BOLD",
          fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
summary <- data.frame(Sheet=c("Protein_0.01","Protein_dr_0.01","Peptides_0.01",
                              "Protein_0.01_0.001","Protein_dr_0.01_0.001","Peptides_0.01_0.001",
                              "Protein","Protein_dr","Peptides","Replication","deCODE","UKB_PPP"),
                      Description=c("Unfiltered proteins, MAF>=0.01","DR-filtered proteins, MAF>=0.01","All peptides, MAF>=0.01",
                                    "Unfiltered proteins, MAF in (0.01,0.001]","DR-filtered proteins, MAF (0.01,0.001]",
                                    "All peptides, MAF in (0.01,0.001]",
                                    "Unfiltered proteins","DR-filtered proteins","All peptides","All/DR-filtered replication",
                                    "deCODE replication","UKB-PPP replication"))
writeDataTable(wb, "Summary", summary, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
for (i in c("Protein_0.01","Protein_dr_0.01","Peptides_0.01",
            "Protein_0.01_0.001","Protein_dr_0.01_0.001","Peptides_0.01_0.001",
            "Protein","Protein_dr","Peptides","Replication","deCODE","UKB_PPP"))
{
    sheetnames <- i
    cat(sheetnames,"\n")
    addWorksheet(wb, sheetnames, zoom=150)
    writeData(wb, sheetnames, sheetnames, xy=c(1,1),
                  headerStyle=createStyle(textDecoration="BOLD",
                                          fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
    body <- get(i)
    writeDataTable(wb, sheetnames, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
    freezePane(wb, sheetnames, firstCol=TRUE, firstActiveRow=3)
}
bStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
hStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
END
}

tables
