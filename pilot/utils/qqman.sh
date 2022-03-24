# 5-4-2021 JHZ

function plotI()
{
  parallel -j1 --env caprion -C' ' '
    export uniprot={1}
    export protein={2}
    gunzip -c bgen/${uniprot}-plink2.gz | \
    cut -f1,2,12 --output-delimiter=" " > ${caprion}/EPCR-PROC/${uniprot}.txt
    R --slave --vanilla --args \
      input_data_path=${caprion}/EPCR-PROC/${uniprot}.txt \
      output_data_rootname=${uniprot}_turboqq \
      plot_title="${protein} (${uniprot})" < ~/cambridge-ceu/turboqq/turboqq.r
    R --slave --vanilla --args \
      input_data_path=${caprion}/EPCR-PROC/${uniprot}.txt \
      output_data_rootname=${uniprot}_turboman \
      custom_peak_annotation_file_path=${caprion}/EPCR-PROC/${uniprot}.annotate \
      reference_file_path=~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda \
      pvalue_sign=8.210181e-12 \
      plot_title="${protein} (${uniprot})" < ~/cambridge-ceu/turboman/turboman.r
    rm ${caprion}/EPCR-PROC/${uniprot}.txt
#   R --no-save < utils/qqman.R
  '
}

function plotII()
{
  parallel -j1 --env caprion -C' ' '
    export uniprot={1}
    export protein={2}
    gunzip -c bgen2/${protein}-plink2.gz | \
    cut -f1,2,12 --output-delimiter=" " > ${caprion}/EPCR-PROC/${protein}.txt
    R --slave --vanilla --args \
      input_data_path=${caprion}/EPCR-PROC/${protein}.txt \
      output_data_rootname=${protein}_turboqq \
      plot_title="${protein} (${uniprot})" < ~/cambridge-ceu/turboqq/turboqq.r
    R --slave --vanilla --args \
      input_data_path=${caprion}/EPCR-PROC/${protein}.txt \
      output_data_rootname=${protein}_turboman \
      custom_peak_annotation_file_path=${caprion}/EPCR-PROC/${protein}.annotate \
      reference_file_path=~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda \
      pvalue_sign=8.210181e-12 \
      plot_title="${protein} (${uniprot})" < ~/cambridge-ceu/turboman/turboman.r
    rm ${caprion}/EPCR-PROC/${protein}.txt
#   R --no-save < utils/qqman.R
  '
}

export caprion=~/rds/projects/Caprion_proteomics/pilot/
(
  echo chromosome position
  sed '1,2d' ${caprtion}/EPCR-PROC/EPCR.sentinels | tr '|' '\t' | cut -f1 | \
  awk '{split($1,a,":");sub(/chr/,"",a[1]);print a[1],a[2]}'
) > ${caprion}/EPCR-PROC/Q9UNN8_invn.annotate
(
  echo chromosome position
  sed '1,2d' ${caprion}/EPCR-PROC/PROC.sentinels | tr '|' '\t' | cut -f1 | \
  awk '{split($1,a,":");sub(/chr/,"",a[1]);print a[1],a[2]}'
) > ${caprion}/EPCR-PROC/P04070_invn.annotate

awk 'NR==1{print $1 "_invn", $2 "_invn"}' qqman.list | plotI
awk 'NR==2{print $1 "_invn", $2 "_invn"}' qqman.list | plotI

(
  echo chromosome position
  cut -d' ' -f2,3 --complement ${caprion}/bgen2/EPCR-PROC/5e-8/caprion-invn.sentinels | \
  grep EPCR_All | \
  sed 's/Chrom/chr/;s/_invn//g;s/chr[0-9]*://' | \
  cut -d' ' -f1,4
) > ${caprion}/EPCR-PROC/EPCR_All_invn.annotate

echo Q9UNN8 EPCR | awk 'NR==1{print $1 "_All_invn", $2 "_All_invn"}' | plotII
