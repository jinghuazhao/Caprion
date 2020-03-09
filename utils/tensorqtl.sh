# 9-3-2020 JHZ

# recycle tensorQTL venv
# ln -s $HOME/tensorqtl/venv
# pip install tensorqtl

module load python/3.6
virtualenv --system-site-package venv
source venv/bin/activate
python -m pip install rpy2
hostname
jupyter notebook --ip=127.0.0.1 --no-browser --port 8087
# ssh -4 -L 8087:127.0.0.1:8087 -fN login-e-10.hpc.cam.ac.uk
firefox http://127.0.0.1:8087/?token=d991ea12ce42b216d3aacd3c573e95280b6cd30d4b4aeeed &

export covariates_file=1.covariates.txt
export prefix=caprion

ln -sf caprion.01.bed 1.bed
awk -v OFS="\t" '{$1="chr" $1};1' caprion.01.bim > 1.bim
ln -sf caprion.01.fam 1.fam
gunzip -c caprion.expression.bed.gz | grep -v chrX | gzip -f > 1.expression.bed.gz

export plink_prefix_path=1
export expression_bed=1.expression.bed.gz

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis

export plink_prefix_path=caprion.01
export expression_bed=caprion.expression.bed.gz

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode trans

R --no-save -q <<END
  library(SparkR)
  sparkR.session()
  df <- read.parquet("caprion.trans_qtl_pairs.parquet")
  head(df)
END
