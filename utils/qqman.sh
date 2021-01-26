# 5-3-2020 JHZ

cat qqman.list | \
parallel -j1 -C' ' '
  echo {1} {2}
  export uniprot={1}
  export protein={2}
  for v in ${uniprot} ${uniprot}_invn
  do
    gunzip -c bgen/${v}-plink2.gz | \
    cut -f1,2,12 --output-delimiter=" " > ${v}.txt

    R --slave --vanilla --args \
      input_data_path=${v}.txt \
      output_data_rootname=${v}_qq \
      plot_title="${v} (${protein})" < turboqq.r

#    R --slave --vanilla --args \
#      input_data_path=${v}.txt \
#      output_data_rootname=${v}_man \
#      custom_peak_annotation_file_path=annotate.txt \
#      reference_file_path=turboman_hg19_reference_data.rda \
#      pvalue_sign=8.210181e-12 \
#      plot_title="${v} (${protein})" < turboman.r

    rm ${v}.txt
  done
'
# Feed to Tryggve
# cat qqman.list | parallel --dry-run -C' ' 'get {1}_invn-plink2.gz'
#
cat qqman.list | parallel -j1 -C' ' 'export uniprot={1};R --no-save < utils/qqman.R'
cat qqman.list | parallel -j1 -C' ' 'export uniprot={1}_invn;R --no-save < utils/qqman.R'
