# 12-12-2019 JHZ

export uniprot=(P12318 P12318 P05362 Q6P179 P0DJI8 P04004 O75347 P51888 P35247)
export protein=(FCGR2A FCGR2A ICAM1 ERAP2 SAA1 VTN TBCA PRELP SFTPD)
export rsid=(rs1801274 rs148396952 rs5498 rs2927608 rs35179000 rs704 rs429358 rs1138545 rs62143206)
export snpid=(1:161479745 1:161569951 19:10395683 5:96252432 11:18290903 17:26694861 19:45411941 9:117835899 19:54326212)
export mi=merged_imputation

function samples()
{
  (
    echo ID_1 ID_2 missing
    echo 0 0 0
    cut -d' ' -f1-3 ${mi}.fam
  ) > ${mi}.sample
  qctool -g ${mi}.bed -s ${mi}.sample -sample-stats -osample ${mi}.sample-stats
  awk 'NR>10 && !/success/' ${mi}.sample-stats | cut -f1,3 > ${mi}.missing
  plink2 --bfile merged_imputation --pca 20 --threads 4 --out merged_imputation
  cut -d' ' -f1,2 SomaLogic.sample | \
  sed '1,2d' | \
  grep -v -f affymetrix.id - > affymetrix.id2
}

function bgen_gen()
{
  for i in `seq 0 8`
  do
    export chr=$(awk -v snpid=${snpid[i]} 'BEGIN{split(snpid,chrpos,":");print chrpos[1]}')
    echo ${rsid[i]} > ${rsid[i]}
    qctool -g impute_${chr}_interval.bgen -s SomaLogic.sample -incl-rsids ${rsid[i]} -incl-samples affymetrix.id -ofiletype bgen_v1.1 -og ${rsid[i]}.bgen
    bgenix -g ${rsid[i]}.bgen -index -clobber
    qctool -g impute_${chr}_interval.bgen -s SomaLogic.sample -incl-rsids ${rsid[i]} -incl-samples affymetrix.id -ofiletype gen -og ${rsid[i]}.gen
    ( head -2 SomaLogic.sample; sed '1,2d' SomaLogic.sample | grep -v "NA NA NA" ) > ${rsid[i]}.sample
  done
}

function assoc_bolt()
{
  bolt \
      --bfile=${mi} \
      --remove=affymetrix.id2 \
      --bgenFile=${rsid[i]}.bgen \
      --sampleFile=${rsid[i]}.sample \
      --phenoFile=SomaLogic.sample \
      --phenoCol=${uniprot[i]} \
      --covarFile=SomaLogic.sample \
      --covarCol=sex \
      --qCovarCol=age \
      --qCovarCol=bmi \
      --qCovarCol=PC{1:20} \
      --lmm \
      --LDscoresUseChip \
      --noMapCheck \
      --numLeaveOutChunks 2 \
      --statsFileBgenSnps=${uniprot[i]}-${rsid[i]}.bgen-stats \
      --statsFile=${uniprot[i]}-${rsid[i]}.stats \
      2>&1 | tee ${uniprot[i]}-${rsid[i]}-bolt.log
}

function assoc_snptest()
{
  snptest \
          -data ${rsid[i]}.bgen ${rsid[i]}.sample -log ${uniprot[i]}-${rsid[i]}-snptest.log -cov_all \
          -filetype bgen \
          -frequentist 1 -hwe -missing_code NA,-999 -use_raw_covariates -use_raw_phenotypes \
          -method score \
          -pheno ${uniprot[i]} -printids \
          -o ${uniprot[i]}-${rsid[i]}.out
  snptest \
          -data ${rsid[i]}.bgen ${rsid[i]}.sample -log ${uniprot[i]}_invn-${rsid[i]}-snptest.log -cov_all \
          -filetype bgen \
          -frequentist 1 -hwe -missing_code NA,-999 -use_raw_covariates -use_raw_phenotypes \
          -method score \
          -pheno ${uniprot[i]}_invn -printids \
          -o ${uniprot[i]}_invn-${rsid[i]}.out
}

for i in `seq 0 8`
do
  assoc_bolt
  assoc_snptest
done
(
  cat *out | head -19 | sed 's/allele//g;s/frequentist_//g' | tail -n 1 | awk -v OFS="\t" '{print "uniprot", "protein", $0}'
  for i in `seq 0 8`
  do
    awk 'NR==20' ${uniprot[i]}-${rsid[i]}.out | \
    awk -v uniprot=${uniprot[i]} -v protein=${protein[i]} -v OFS="\t" '{print uniprot, protein, $0}'
  done
) > affymetrix.tsv
(
  cat *out | head -19 | sed 's/allele//g;s/frequentist_//g' | tail -n 1 | awk -v OFS="\t" '{print "uniprot", "protein", $0}'
  for i in `seq 0 8`
  do
    awk 'NR==20' ${uniprot[i]}_invn-${rsid[i]}.out | \
    awk -v uniprot=${uniprot[i]} -v protein=${protein[i]}_invn -v OFS="\t" '{print uniprot, protein, $0}'
  done
) > affymetrix_invn.tsv

R --no-save -q < ps.R

(
  head -1 SOMALOGIC_Master_Table_160410_1129info.tsv
  cut -f1 affymetrix.tsv | sed '1d' | grep -w -f - SOMALOGIC_Master_Table_160410_1129info.tsv
) > SomaLogic.info

join -j2 <(cut -f2,10 SomaLogic.info | sed '1d'| sort | uniq) \
         <(cut -f10 SomaLogic.info | sed '1d'| sort | uniq | grep -f - -w glist-hg19 | cut -d' ' -f1,4) | \
cut -d' ' -f2,3 | \
parallel -j1 -C' ' 'cd SomaLogic;unzip {1}.zip {1}/{1}_chrom_{2}_meta_final_v1.tsv.gz;cd -'

(
  echo "SOMAS_ID_round1 chr pos rsid protein VARIANT_ID chromosome position Allele1 Allele2 Effect StdErr logP"
  join -j2 <(cut -f2,10 SomaLogic.info | sed '1d'| sort | uniq) \
           <(cut -f10 SomaLogic.info | sed '1d'| sort | uniq | grep -f - -w glist-hg19 | cut -d' ' -f1,4) | \
  join - <(cut -f2,3,6 affymetrix.tsv | sed '1d' | sort -k1,1) | cut -d' ' -f1-5 | \
  parallel -j1 -C' ' 'echo {2} {3} {5} {4} {1} $(zgrep -w {5} SomaLogic/{2}/{2}_chrom_{3}_meta_final_v1.tsv.gz)'
) > SomaLogic.box

R --no-save -q <<END
  s <- read.delim("SomaLogic.info",as.is=T)
  v <- c("chromosome_number","start_position","end_position","ensembl_gene_id","external_gene_name","UniProt","Target")
  s <- within(s[v],{chromosome_number <- paste0("chr",chromosome_number)})
  names(s) <- c("chrom","chromStart","chromEnd","ENSEMBL","HGNC","uniprot","Protein")
  write.table(s,file="SomaLogic.gene",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
END

(
  echo $(head -1 SomaLogic.gene | sed 's/chrom/geneChr/;s/chromStart/geneStart/;s/chromEnd/geneEnd/') \
       $(head -1 affymetrix.results | sed 's/chrom/snpChr/;s/chromStart/snpStart/;s/chromEnd/snpEnd/') | \
  sed 's/ /\t/g'
  bedtools intersect -a SomaLogic.gene -b affymetrix.results -wo
) | \
cut -f1-20 > SomaLogic.tsv

(
  echo chrom chromStart chromEnd gene
  cut -f2 affymetrix.tsv | \
  sed '1d' | \
  grep -w -f - glist-hg19 | \
  awk -v OFS="\t" '{print "chr" $1,$2,$3,$4}'
) | \
bedtools intersect -a - -b affymetrix.results -wo 

# ERAP2 peptides and rs2927608
export rsid=rs2927608
for peptides in $(head -1 ERAP2.sample | sed 's/ /\n/g' | awk 'NR>26')
do
  snptest \
          -data ${rsid}.bgen ERAP2.sample -log ERAP2-${peptides}-${rsid}-snptest.log -cov_all \
          -filetype bgen \
          -frequentist 1 -hwe -missing_code NA,-999 -use_raw_covariates -use_raw_phenotypes \
          -method score \
          -pheno ${peptides} -printids \
          -o ERAP2-${peptides}-${rsid}-snptest.out
done

(
  awk 'NR==19' ERAP2-*-${rsid}-snptest.out | head -1 | awk '{print "peptides_Isotope.Group", $0}'
  for peptides in $(head -1 ERAP2.sample | sed 's/ /\n/g' | awk 'NR>26')
  do
     awk -v peptides=${peptides} 'NR==20 {print peptides, $0}' ERAP2-${peptides}-${rsid}-snptest.out
  done
) | \
sed 's/ /\t/g;s/frequentist_//g' > ERAP2.tsv
