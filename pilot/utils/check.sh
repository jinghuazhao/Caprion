# 3-3-2020

join <(cut -d',' -f1,2 protein_list.csv| awk -vFS="," '{print $2,$1}' | sort -k1,1) \
     <(sed 's/ /\n/g' ../SWATH-MS/swath-ms.uniprot | sort -k1,1) > caprion-swath-ms.overlap
