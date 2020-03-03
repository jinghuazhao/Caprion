# 3-3-2020 JHZ

gunzip -c plink2/P04004-plink2.gz | cut -f1,2,12 --output-delimiter=' ' > VTN.txt

R --slave --vanilla --args \
  input_data_path=VTN.txt \
  output_data_rootname=VTN_qq \
  plot_title="VTN example" < turboqq.r

R --slave --vanilla --args \
  input_data_path=VTN.txt \
  output_data_rootname=VTN_man \
  reference_file_path=turboman_hg19_reference_data.rda \
  pvalue_sign=5e-8 \
  plot_title="VTN example" < turboman.r

