#/usr/bin/bash

export TMPDIR=$HPC_WORK/work
export tag=_nold

function nogrouping()
{
  for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do
  echo $p
  export p=${p}
  (
    mergeBed -i sentinels/${p}_nold.p -d 1000000 -c 13 -o min | \
    awk -v OFS="\t" -v prot=${p} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "prot"
      print $0, prot
    }'
  ) > work/${p}.merged
  (
    cut -f1-4,13 sentinels/${p}_nold.p | \
    bedtools intersect -a work/${p}.merged -b - -wa -wb | \
    awk '$4==$10' | \
    cut -f1-6,8-10 | \
    awk -v OFS="\t" '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "prot", "MarkerName", "CHR", "POS", "SNP", "P_check"
      $5=$5 OFS $6 ":" $7
      gsub(/chr/,"",$6)
      print
    }'
  ) | uniq > work/${p}.sentinels
  done
  (
    cat work/*sentinels | head -1
    for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do awk 'NR>1' work/${p}.sentinels; done
  ) > caprion.merge
  cut -f5 caprion.merge | sed '1d' | sort | uniq > caprion.merge.prot
  R --no-save -q <<\ \ END
    merge <- within(read.delim("caprion.merge",as.is=TRUE),{Chrom <- as.numeric(gsub("chr","",Chrom))})
    m <- subset(merge,SNP!=".")
    ord <- with(m,order(Chrom,End))
    write.table(m[ord,c(1:6,9)],file="caprion-invn.sentinels",row.names=FALSE,quote=FALSE)
    library(TwoSampleMR)
    options(width=200)
    write.table(round(ld_matrix(c("rs867186",subset(m,Chrom==20)[["SNP"]]),with_alleles=FALSE),digits=3),file="caprion-20ld.tsv",quote=FALSE,sep="\t")
  END
  # Not in the LD reference panel
  # rs545281213, rs190767097, rs117612661, rs33970207, rs139746733
  cut -d' ' -f5 caprion-invn.sentinels | sed '1d' | sort | uniq > caprion-invn.sentinels.prot
  gunzip -c hgTables.gz | awk 'length($1)<=5' | grep -f caprion-invn.sentinels.prot -
  # chr20:33764554_A_G rs867186
  cd ${caprion}/bgen2
  (
    gunzip -c *.gz | \
    head -1
    zgrep -w rs867186 *gz
  ) | sed 's/_invn-plink2.gz://g' > ${caprion}/rs867186.txt
  cd -
}

function grouping()
{
  for p in $(ls sentinels/*${group}${tag}.p | grep ${group} | sed 's|sentinels/||g;s|'"-${group}$tag"'.p||g'); do
  echo $p
  export p=${p}
  (
    mergeBed -i sentinels/${p}-${group}_nold.p -d 1000000 -c 13 -o min | \
    awk -v OFS="\t" -v prot=${p} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "prot"
      print $0, prot
    }'
  ) > work/${p}-${group}.merged
  (
    cut -f1-4,13 sentinels/${p}-${group}_nold.p | \
    bedtools intersect -a work/${p}-${group}.merged -b - -wa -wb | \
    awk '$4==$10' | \
    cut -f1-6,8-10 | \
    awk -v OFS="\t" '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "prot", "MarkerName", "CHR", "POS", "SNP", "P_check"
      $5=$5 OFS $6 ":" $7
      gsub(/chr/,"",$6)
      print
    }'
  ) | uniq > work/${p}-${group}.sentinels
  done
  (
    cat work/*${group}.sentinels | head -1
    for p in $(ls sentinels/*${group}${tag}.p | sed 's|sentinels/||g;s|'"-${group}$tag"'.p||g'); do awk 'NR>1' work/${p}-${group}.sentinels; done
  ) > caprion-${group}.merge
  cut -f5 caprion-${group}.merge | sed '1d' | sort | uniq > caprion-${group}.merge.prot
  R --no-save -q <<\ \ END
    group <- Sys.getenv("group")
    merge <- within(read.delim(paste0("caprion-",group,".merge"),as.is=TRUE),{Chrom <- as.numeric(gsub("chr","",Chrom))})
    m <- subset(merge,SNP!=".")
    ord <- with(m,order(Chrom,End))
    write.table(m[ord,c(1:6,9)],file=paste0("caprion-",group,"-invn.sentinels"),row.names=FALSE,quote=FALSE)
    library(TwoSampleMR)
    options(width=200)
    write.table(round(ld_matrix(c("rs867186",subset(m,Chrom==20)[["SNP"]]),with_alleles=FALSE),digits=3),
                file=paste0("caprion-",group,"-20ld.tsv"),quote=FALSE,sep="\t")
  END
  # Not in the LD reference panel
  # group1 rs146082057, rs544527276, rs145965345, rs527241169, rs397865498
  # group2 rs142132444
  cut -d' ' -f5 caprion-${group}-invn.sentinels | sed '1d' | sort | uniq > caprion-${group}-invn.sentinels.prot
  gunzip -c hgTables.gz | awk 'length($1)<=5' | grep -f caprion-${group}-invn.sentinels.prot -
  # chr20:33764554_A_G rs867186
  cd ${caprion}/bgen2
  (
    gunzip -c *-${group}-plink2.gz | \
    head -1
    zgrep -w rs867186 *${group}-plink2.gz
  ) | sed 's/_invn-${group}-plink2.gz://g' > ${caprion}/rs867186-${group}.txt
  cd -
}

export group=group1
grouping
export group=group2
grouping
