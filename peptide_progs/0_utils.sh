#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export signals=${analysis}/work/caprion${suffix}.signals

# All signals
cat <(cat ${analysis}/peptide/*/*signals | head -1 | paste <(echo protein) -) \
    <(ls ${analysis}/peptide/*/*signals | xargs -l basename -s .signals | grep -e PROC -e EPCR -e INHBE -v | \
      parallel -C' ' 'export prot={};awk -vOFS="\t" "NR>1{print ENVIRON[\"prot\"],\$0}" ${analysis}/peptide/{}/{}.signals') | \
      awk -vOFS="\t" '{print $1,$2,$6,$8}' > ${analysis}/reports/peptide.signals
# cis/trans
cat <(cat ${analysis}/peptide/*/*cis.vs.trans | head -1) \
    <(ls ${analysis}/peptide/*/*cis.vs.trans | xargs -l basename -s .cis.vs.trans | grep -e PROC -e EPCR -e INHBE -v | \
      parallel -C' ' 'awk "NR>1" ${analysis}/peptide/{}/{}.cis.vs.trans') \
    > ${analysis}/reports/peptide.cis.vs.trans

# Proteins
sed '1d' ${analysis}/reports/peptide.signals | \
cut -f1 | \
sort | \
grep -e PROC -e EPCR -e INHBE -v | \
uniq -c | \
wc -l

# Left-over
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

# Protein-peptides
sed '1d' ${analysis}/reports/peptide.signals | \
cut -f1,2 | \
sort | \
uniq | \
wc -l

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

# R counterpart
# Note also fasta_sequences <- readDNAStringSet(fasta_file_path, format = "fasta")
R --no-save <<END
library(Biostrings)
fasta_file_path <- 'P04217.fasta'
search_442688365 <- 'TDGEGALSEPSATVTIEELAAPPPPVLMHHGESSQVLHPGNK'
fasta_sequences <- readAAStringSet(fasta_file_path, format = "fasta")
first_sequence <- fasta_sequences[[1]]
cat("Sequence:", toString(first_sequence), "\n")
match_position <- regexpr(search_442688365, first_sequence)
match_position
END
