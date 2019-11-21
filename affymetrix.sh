# 20-11-2019 JHZ

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
  (
    cut -d ' ' -f1-2,4-36 SomaLogic.sample | \
    sed 's/ID_1/FID/g;s/ID_2/IID/g;2d'| grep -v "NA NA NA"
  ) > SomaLogic.covar.pheno
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
      --phenoFile=SomaLogic.covar.pheno \
      --phenoCol=${uniprot[i]} \
      --covarFile=SomaLogic.covar.pheno \
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
