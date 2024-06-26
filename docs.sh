#!/usr/bin/bash

function setup()
{
  module load python/3.8
  source ~/rds/public_databases/software/py38/bin/activate
  for f in .gitignore README.md autoencoder.md caprion.ini caprion.ipynb caprion.R docs.sh gwas2.md mkdocs.yml utils
  do cp -r -p ../pilot/${f} pilot; done
}

if [ "$(uname -n | sed 's/-[0-9]*$//')" == "login-q" ]; then
   module load ceuadmin/libssh/0.10.6-icelake
   module load ceuadmin/openssh/9.7p1-icelake
fi

setup
mkdocs build
mkdocs gh-deploy

git add workflow
git commit -m "Snakemake profiles"
git add .gitignore
git commit -m ".gitignore"
git add README.md pilot/README.md
git commit -m "README"
git add pilot/qrcode.svg pilot/qrcode.png
git commit -m "QR code"
git add peptide_progs/ 0_compare.sb 0_utils.sb
git add 1_pca_projection.sh 2_ggm.R 3_wgcna.R 4_pca_clustering.R 5_pgwas.* 6_meta_analysis.* 7_merge* 8_hla.* 9_lookup.* ToDo.sh
git commit -m "Analysis"
git add pilot/caprion.ini pilot/caprion.R pilot/caprion.ipynb
git commit -m "Primary Code"
git add pilot/gwas2.md
git commit -m "gwas2"
git add pilot/autoencoder.md
git commit -m "autoencoder"
git add docs.sh pilot/.gitignore pilot/utils pilot/docs.sh pilot/Notes.md
git commit -m "utils"
git add mkdocs.yml pilot/mkdocs.yml
git commit -m "mkdocs.yml"
git push
# git remote set-url origin https://github.com/jinghuazhao/Caprion.git
