#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export signals=${analysis}/work/caprion${suffix}.signals

module load ceuadmin/R

# All signals
cat <(cat ${analysis}/peptide/*/*signals | head -1 | paste <(echo protein) -) \
    <(ls ${analysis}/peptide/*/*signals | xargs -l basename -s .signals | \
      parallel -C' ' 'export prot={};awk -vOFS="\t" "NR>1{print ENVIRON[\"prot\"],\$0}" ${analysis}/peptide/{}/{}.signals') | \
      awk -vOFS="\t" '{print $1,$2,$6,$8}' > ${analysis}/reports/peptide.signals
# cis/trans
cat <(cat ${analysis}/peptide/*/*cis.vs.trans | head -1) \
    <(ls ${analysis}/peptide/*/*cis.vs.trans | xargs -l basename -s .cis.vs.trans | \
      parallel -C' ' 'awk "NR>1" ${analysis}/peptide/{}/{}.cis.vs.trans') \
    > ${analysis}/reports/peptide.cis.vs.trans

Rscript -e '
  suppressMessages(library(dplyr))
  suppressMessages(library(gap))
  analysis <- getwd()
  cis.vs.trans <- read.csv(file.path(analysis,"reports","peptide.cis.vs.trans"))
  png(file.path(analysis,"reports",paste0("peptide",".pqtl2d.png")),width=12,height=10,unit="in",res=300)
  r <- qtl2dplot(cis.vs.trans,chrlen=gap::hg19,snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                 gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                 TSS=TRUE,cis="cis",plot=TRUE,cex.labels=0.6,cex.points=0.6,
                 xlab="pQTL position",ylab="Gene position")
  dev.off()
  r <- qtl2dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                   snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                   gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                   TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                   xlab="pQTL position",ylab="Gene position")
  htmlwidgets::saveWidget(r,file=file.path(analysis,"reports",paste0("peptide",".pqtl2dplotly.html")))
  r <- qtl3dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                   snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                   gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                   TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                   xlab="pQTL position",ylab="Gene position")
  htmlwidgets::saveWidget(r,file=file.path(analysis,"reports",paste0("peptide",".pqtl3dplotly.html")))
  snplist_01 <- scan("~/Caprion/analysis/bgen/caprion-0.01.snplist",what="")
  cvt_01 <- read.csv("~/Caprion/analysis/reports/peptide.cis.vs.trans") %>%
            filter(SNP %in% snplist_01)
  write.csv(cvt_01,file="~/Caprion/analysis/reports/peptide-01.cis.vs.trans",row.names=FALSE,quote=FALSE)
'

# Proteins
sed '1d' ${analysis}/reports/peptide.signals | \
cut -f1 | \
sort | \
uniq -c | \
wc -l

Rscript -e '
  signals <- read.delim("~/Caprion/analysis/reports/peptide.signals")
  dim(signals)
  tbl <- with(signals,table(protein))
  png("~/Caprion/analysis/reports/peptide.png",width=6,height=4,res=300,units="in")
  hist(tbl,main="Number of proteins by signal",xlab="No. signals", ylab="No. proteins")
  dev.off()
'

# Left-over
## proteins
(
  export n_with_signals=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | wc -l)
  for i in $(seq 300) # $(seq ${n_with_signals})
  do
    export signal_index=${i}
    export protein=$(awk 'NR>1{print $1}' ${signals} | sort -k1,1 | uniq | awk 'NR==ENVIRON["signal_index"]')
    export root=${analysis}/peptide/${protein}
    export pheno=${root}/${protein}.pheno
    export N=$(awk 'NR==1{print NF-2}' ${pheno})
    export all_peptides=$(head -1 ${pheno} | cut -f1,2 --complement)
    export pqtl_peptides=$(sed '1d' ${root}/${protein}.signals | cut -f1 | sort -k1,1n | uniq)
    export array=$(grep -n -f <(echo ${pqtl_peptides} | tr ' ' '\n') <(echo ${all_peptides} | tr ' ' '\n') | cut -d':' -f1 | tr '\n' ',' | sed 's/.$//')
    export dir=${root}/qqmanhattanlz
    echo ${signal_index}, ${protein}
  done
) | \
grep -f <(sed '1d' ${analysis}/reports/peptide.signals | cut -f1 | sort | uniq) -v - | cut -d',' -f1 > ${analysis}/left-over

## peptides
## CO3, ITIH2 dosage>>genotype since dosage.png is complete but its raw files may be missing as well
ls *dosage.png | sed 's/-dosage.png//'  | grep -v -f <(ls *genotype.png | sed 's/-genotype.png//'| sort -k1,1) | sed 's/-/\t/g' | cut -f2 | uniq
ls *.dat | sed 's/.dat//'  | grep -v -f <(ls *.raw | sed 's/.raw//'| sort -k1,1) | sed 's/-/\t/g' | cut -f3 | uniq

# Protein-peptides
sed '1d' ${analysis}/reports/peptide.signals | \
cut -f1,2 | \
sort | \
uniq | \
wc -l

# peptide-pqtl
cd CO3
cut -f1,7 CO3.signals | sed '1d;s/\t/_/;s/X:[0-9]+_[A-Z]+[A-Z]+/X:[0-9]+/' | sort | \
join -v1 - <(ls qqmanhattanlz/*svg | xargs -l basename -s .svg | sed 's/chr23_/X:/' | sort) | cut -d'_' -f1 | sort | uniq | tr '\n' ','

cd ${analysis}/peptide
for prot in APOB EPCR PROC ERAP2 ITIH2
do
  echo $prot
  echo invnormal results
  cut -f2 $prot/$prot.signals  | sort -n | uniq -c | awk 'NR==1{print "Chr|pQTL"}{print $2"|"$1}'
  echo scaled results
  cut -f2 01-15-2024-keep/$prot/$prot.signals  | sort -n | uniq -c | awk 'NR>1{print $2"|"$1}'
done
cd -

# discrepancy in signals vs cis/trans classification
# GP1BA missing
export f1=${analysis}/reports/peptide.signals
export f2=${analysis}/reports/peptide.cis.vs.trans
awk -vFS=, 'NR>1{print $3"-"$5"-"$2}' ${f2} | sort | paste - <(sed '1d' ${f1} | cut -f1,2,4 | tr '\t' '-' | sort) | awk -vFS="\t" '$1!=$2'

function bgz()
{
  module load samtools/1.13/gcc/zwxn7ug3
  if [ ! -d ${analysis}/METAL${suffix}/gz ]; then mkdir ${analysis}/METAL${suffix}/gz; fi
  ls ${analysis}/METAL${suffix}/*-1.tbl.gz | \
  xargs -l basename -s -1.tbl.gz | \
  parallel -j10 -C' ' '
  echo {}
  (
    echo chromsome position log_pvalue beta se
    gunzip -c ${analysis}/METAL${suffix}/{}-1.tbl.gz | \
    awk "{print \$1,\$2,-\$12,\$10,\$11}"
  ) | \
  bgzip -f > ${analysis}/METAL${suffix}/gz/{}.txt.gz
  '
}

function chrX()
{
cut -d' ' -f1 ${analysis}/output/caprion${suffix}-1.id | grep -f - ${analysis}/output/chrX.idlist \
   > ${analysis}/output/chrX${suffix}-1.id
cut -d' ' -f1 ${analysis}/output/caprion${suffix}-2.id | grep -f - ${analysis}/output/chrX.idlist \
   > ${analysis}/output/chrX${suffix}-2.id
cut -d' ' -f1 ${analysis}/output/caprion${suffix}-3.id | grep -f - ${analysis}/output/chrX.idlist \
   > ${analysis}/output/chrX${suffix}-3.id
}

function pav()
{
  for protein in A1BG APOB EPCR ERAP2 PROC PON3
  do
  export protein=${protein}
  export root=${analysis}/peptide/${protein}
  cd ${root}/METAL/vep
  (
    grep -v '##' *.tab | head -1
    grep -v '#' *.tab
  ) > ${root}/vep.txt
  cd - > /dev/null
  Rscript -e '
   analysis <- Sys.getenv("analysis")
   root <- Sys.getenv("root")
   protein <- Sys.getenv("protein")
   vep <- read.delim(file.path(root,"vep.txt"),check.names=FALSE)
   knitr::kable(with(vep,table(Consequence)),caption=paste(protein,"annotation"))
  '
  done
  (
    head -1 ${analysis}/peptide/*/vep.txt | awk "/Consequence/" | head -1
    grep -v '#' ${analysis}/peptide/*/vep.txt | sed "s|${analysis}/peptide/||;s|vep.txt:||;s|.tab:|,|"
  ) > ${analysis}/work/peptide_vep.txt
  Rscript -e '
   analysis <- Sys.getenv("analysis")
   vep <- read.delim(file.path(analysis,"work","peptide_vep.txt"),check.names=FALSE)
   knitr::kable(with(vep,table(Consequence)),caption="Peptide annotation based on six proteins")
  '
}

function meta()
{
  Rscript -e '
     options(width=200)
     suppressMessages(library(Biobase))
     suppressMessages(library(dplyr))
     suppressMessages(library(tidyr))
     suppressMessages(library(stringr))
     suppressMessages(library(Biostrings))
     analysis <- "~/Caprion/analysis"
     load("~/Caprion/pilot/ZYQ.rda")
     isotopes <- vector()
   # for (p in c("A1BG","APOB","EPCR","ERAP2","PROC"))
     p <- "PROC"
     protein_names <- paste0(p,"_HUMAN")
     benchmarks <- grep(paste(protein_names,collapse="|"),names(uniprot_db))
     {
        isotope <- read.table(file.path(analysis,"peptide",p,paste0(p,".pheno")),
                              check.names=FALSE,header=TRUE, nrows=1) %>%
                   dplyr::select(-FID,-IID) %>%
                   names()
        isotopes <- c(isotopes,isotope)
     }
     peptide_data <- mapping_ZYQ %>%
                     dplyr::filter(Isotope.Group.ID %in% isotopes) %>%
                     dplyr::select(Isotope.Group.ID,Modified.Peptide.Sequence,Monoisotopic.m.z,Max.Isotope.Time.Centroid,Charge)
     fasta_file <- file.path(analysis,"crux","uniprot_sprot.fasta")
     uniprot_db <- readAAStringSet(fasta_file)
     unique_peptides <- peptide_data %>%
       group_by(Isotope.Group.ID) %>%
       summarise(
         Modified.Peptide.Sequence = dplyr::first(Modified.Peptide.Sequence),
         Monoisotopic.m.z = mean(Monoisotopic.m.z),
         Charge = dplyr::first(Charge),
         Max.Isotope.Time.Centroid = mean(Max.Isotope.Time.Centroid)
       ) %>%
       ungroup() %>%
       distinct(Modified.Peptide.Sequence, .keep_all = TRUE)
     identify_peptide <- function(peptide, uniprot_db, protein_names) {
       # Remove modifications for exact matching
       peptide_clean <- gsub("\\[.*?\\]", "", peptide)
       hits <- vmatchPattern(peptide_clean, uniprot_db)
       if (length(hits) > 0) {
         matching_proteins <- names(hits)
         matching_protein_names <- protein_names[matching_proteins]
         return(paste(matching_protein_names, collapse = "; "))
       } else {
         return(NA)
       }
     }
     match_peptides <- unique_peptides
     match_peptides$Matched_Proteins <- sapply(unique_peptides$Modified.Peptide.Sequence, function(peptide) {
       identify_peptide(peptide, uniprot_db, protein_names)
     })
     print(match_peptides)
  '
}
Rscript -e '
  peptides <- function(code)
  {
    suppressMessages(library(dplyr))
    cat(code,"\n")
    load(paste0("~/Caprion/pilot/",code,".rda"))
    peptides_per_protein <- get(paste0("mapping_",code)) %>%
                            filter(Protein!="-") %>%
                            group_by(Protein) %>%
                            summarize(N=n(),
                                      group=cut(N, breaks = c(0, 15, 30, 50, Inf),
                                                   labels = c("1-15", "16-30", "31-50", "51+"),
                                                   right = FALSE))
    print(dim(peptides_per_protein))
    attach(peptides_per_protein)
    print(table(N))
    table_group <- table(group)
    table_sum <- table(group)|>sum()
    detach(peptides_per_protein)
    wtable <- with(peptides_per_protein, tapply(N,group,sum))
    list(table_group=c(code,table_group),table_sum=c(code,table_sum),wtable=c(code,wtable))
  }
  ZWK <- peptides("ZWK")
  ZYQ <- peptides("ZYQ")
  UDP <- peptides("UDP")
  UHZ <- peptides("UHZ")
  knitr::kable(rbind(ZWK$table_group,ZYQ$table_group,UDP$table_group,UHZ$table_group),caption="Frequencies of categories")
  knitr::kable(rbind(ZWK$wtable,ZYQ$wtable,UDP$wtable,UHZ$wtable),caption="Total number of peptides")
 '

module load ceuadmin/htslib
export PERL5LIB=

function sumstats()
{
  export prot=${1}
  cd ~/Caprion/analysis
# uncomment below for a fresh run
# ls peptide/${prot}/*gz | parallel -j10 -C' ' 'tabix -f -S1 -s1 -b3 -e3 {}'
  if [ ! -d peptide/${prot}/sumstats ]; then mkdir peptide/${prot}/sumstats; fi
  parallel -j1 -C' ' '
     export isotope_batch=$(basename -s .fastGWA.gz {2})
     echo ${isotope_batch}
     export out=peptide/${prot}/sumstats/${isotope_batch}-{1}.txt
     gunzip -c {2} | \
     head -1 > ${out}
     tabix {2} {1} >> ${out}
  '  ::: $(grep -w ${prot} work/caprion_dr.cis.vs.trans | \
           awk -v FS=, '{print $8,$9}' | \
           awk -vM=1e6 '{print $1":"$2-M"-"$2+M}') ::: $(ls peptide/${prot}/*gz)
  ls peptide/${prot}/sumstats/*.txt | \
  parallel -C' ' 'if [ "$(wc -l < {})" == 1 ]; then rm {}; fi'
  (
    gunzip -c peptide/${prot}/*gz | \
    head -1 | \
    awk -vOFS="\t" '{print "source",$0}'
    (
      parallel -j1 -C' ' '
         export isotope_batch=$(basename -s .fastGWA.gz {2})
         tabix {2} {1} | \
         awk -vOFS="\t" -vsource=${isotope_batch}:{1} "{print source,\$0}"
      '  ::: $(grep -w ${prot} work/caprion_dr.cis.vs.trans | \
               awk -v FS=, '{print $8,$9}' | \
               awk -vM=0 '{print $1":"$2-M"-"$2+M}') ::: $(ls peptide/${prot}/*gz)
    ) | \
    sort -k1,1
  ) > peptide/${prot}/${prot}.sumstats
  cd -
}

sumstats PROC
sumstats ERAP2
sumstats EPCR
sumstats A1BG
sumstats APOB
