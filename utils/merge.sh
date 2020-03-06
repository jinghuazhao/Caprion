# 6-3-2020 JHZ

export TMPDIR=$HPC_WORK/work
export tag=_nold

for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do

echo $p
export p=${p}
(
  mergeBed -i sentinels/${p}_nold.p -d 1000000 -c 13 -o min | \
  awk -v OFS="\t" -v prot=${p} '
  {
    if(NR==1) print "Chrom", "Start", "End", "logP", "prot"
    print $0, prot
  }'
) > work/${p}.merged
(
  cut -f1-4,13 sentinels/${p}_nold.p | \
  bedtools intersect -a work/${p}.merged -b - -wa -wb | \
  awk '$4==$10' | \
  cut -f1-5,9,10 | \
  awk -v OFS="\t" '
  !(/CCL25/&&/chr19:49206145_C_G/){
    if(NR==1) print "Chrom", "Start", "End", "log10p", "prot", "MarkerName", "log10p_check", "CHR", "POS", "SNP", "BP"
    CHR=substr($1,4)
    split($6,noalleles,"_")
    split(noalleles[1],chrpos,":")
    POS=chrpos[2]
    SNP=$6
    BP=chrpos[2]
    print $0,CHR,POS,SNP,BP
  }'
) | uniq > work/${p}.sentinels

done

(
  cat work/*sentinels | head -1
  for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do awk 'NR>1' work/${p}.sentinels; done
) > caprion.merge

R --no-save -q <<END
  merge <- read.delim("caprion.merge",as.is=TRUE)
  m <- subset(merge,SNP!=".")
  write.table(m[,1:6],file="caprion-invn.sentinels",row.names=FALSE,quote=FALSE)
END

cut -d' ' -f5 caprion-invn.sentinels | sed '1d' | sort | uniq > caprion-invn.sentinels.prot
gunzip -c hgTables.gz | awk 'length($1)<=5' | grep -f caprion.sentinels.prot -
