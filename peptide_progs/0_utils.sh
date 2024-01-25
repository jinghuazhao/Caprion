#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export signals=${analysis}/work/caprion${suffix}.signals

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

# discrepancy in signals vs cis/trans classification
# GP1BA missing
export f1=${analysis}/reports/peptide.signals
export f2=${analysis}/reports/peptide.cis.vs.trans
awk -vFS=, 'NR>1{print $3"-"$5"-"$2}' ${f2} | sort | paste - <(sed '1d' ${f1} | cut -f1,2,4 | tr '\t' '-' | sort) | awk -vFS="\t" '$1!=$2'

# A1BG, APOB
# https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28protein+sequence+A1BG_HUMAN%29
# https://rest.uniprot.org/uniprotkb/P04217.fasta
# https://rest.uniprot.org/uniprotkb/P04114.txt

source ~/rds/public_databases/software/py38/bin/activate
python <<END
from Bio import SwissProt
with open('P04217.txt') as file:
    for record in SwissProt.parse(file):
        print(f"ID: {record.entry_name}, Description: {record.description}")
        for feature in record.features:
            print(feature)
        print("=" * 30)
        print("Sequence:")
        print(record.sequence)

from Bio import SeqIO
fasta_file_path = 'P04217.fasta'
sequences = SeqIO.to_dict(SeqIO.parse(fasta_file_path, 'fasta'))
for sequence_id, sequence_record in sequences.items():
    print(f"Sequence ID: {sequence_id}")
    print(f"Sequence Description: {sequence_record.description}")
    print(f"Sequence: {sequence_record.seq}")
    print("=" * 30)

def match_sequence(fasta_file, search_string):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    for sequence_id, sequence_record in sequences.items():
        if search_string in sequence_record.seq:
            print(f"Match found in sequence with ID: {sequence_id}")
            print(f"Sequence Description: {sequence_record.description}")
            print(f"Matched Substring: {search_string}")
            print("=" * 30)
            return True
    print(f"No match found for substring: {search_string}")
    return False

def find_matching_position(sequence, search_string):
    position = sequence.find(search_string)
    return position

fasta_file = 'P04217.fasta'
search_442688365 = 'TDGEGALSEPSATVTIEELAAPPPPVLMHHGESSQVLHPGNK'
match_sequence(fasta_file, search_442688365)

sequence = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
id = 'sp|P04217|A1BG_HUMAN'
amino_acid_sequence = sequence[id].seq
match_position = find_matching_position(amino_acid_sequence, search_442688365)
print(f"Match found at position: {match_position}")
END

#' Peptide mappings to given protein sequence, via
#' 1. accession from panel file, pQTLdata::caprion
#' 2. protein sequence from UniProt
#' 2. modified peptide squence from mapping_ZWK
#' 3. mapping of peptide sequences to protein sequence
#'
#' Note protein name could be with/without _HUMAN suffix & readDNAStringSet
#' It remains to implement a sequence plot with signals

R --no-save <<END
peptideMapping <- function(protein,batch="ZWK",mm=5)
{
  library(Biostrings)
  library(dplyr)
  accession <- subset(pQTLdata::caprion,grepl(protein,Protein))[["Accession"]]
  load(paste0("~/Caprion/pilot/",batch,".rda"))
  mapping <- subset(get(paste0("mapping_",batch)),grepl(protein,Protein))[1:6] %>%
             setNames(c("Isotope.Group.ID","Protein","Modified.Peptide.Sequence", "Monoisotopic.mz", "Max.Isotope.Time.Centroid", "Charge"))
  mps <- mapping %>%
            group_by(Modified.Peptide.Sequence) %>%
            reframe(n_isotope=dplyr::n(),
                   isotope=paste(Isotope.Group.ID,collapse=";"),
                   mz=paste(Monoisotopic.mz,collapse=";"),
                   mtime=paste(Max.Isotope.Time.Centroid,collapse=";"),
                   charge=paste(Charge,collapse=";"))
  fasta_file_path <- paste0("https://rest.uniprot.org/uniprotkb/",accession,".fasta")
  fasta_sequences <- readAAStringSet(fasta_file_path, format = "fasta")
  sequence <- fasta_sequences[[1]]
  cat("Sequence:", toString(sequence), "\n")
  mp <- sapply(setdiff(mapping[["Modified.Peptide.Sequence"]],"-"),
        function(x) {
          y <- gsub("\\[\\d+\\.\\d+\\]", "X", x);
          mp <- matchPattern(y,sequence,with.indels=TRUE,max.mismatch=mm)
          as.data.frame(mp)
        })
  invisible(list(accession=accession,sequence=sequence,mapping=mapping,mps=mps,positions=t(mp)))
}

options(width=200)
A1BG <- peptideMapping("A1BG")
sink("A1BG")
knitr::kable(A1BG$mps[1:3],"markdown")
sink()
unlink("A1BG")
ITIH2 <- peptideMapping("ITIH2")
subset(ITIH2$mapping,Isotope.Group.ID==442581854)
knitr::kable(subset(ITIH2$mapping, rownames(ITIH2$mapping) >=13480 & rownames(ITIH2$mapping) <13492)[c(1,3,4,5,6)],row.names=FALSE)
APOB <- peptideMapping("APOB")
ERAP2 <- peptideMapping("ERAP2",mm=1)

peptideLines <- function(peptides, main="Peptides")
{
  segment_data <- as.data.frame(peptides) %>%
                  mutate(ID = 1:nrow(peptides))
  xmin <- min(unlist(segment_data[["start"]]))
  xmax <- max(unlist(segment_data[["end"]]))
  opar <- par()
  par(mfrow=c(2,1))
  plot(1:10,type="n",axes=FALSE,xlab="",ylab="")
  plot(segment_data$start,segment_data$ID,type="n",ylab="Segment",xlab="Position",xlim=c(xmin,xmax))
  for (i in 1:nrow(segment_data)) {
    lines(c(segment_data$start[i], segment_data$end[i]), c(segment_data$ID[i], segment_data$ID[i]), col = segment_data$ID[[i]], lwd = 2)
  }
  par(opar)
}

peptideLines(ERAP2$positions)

END