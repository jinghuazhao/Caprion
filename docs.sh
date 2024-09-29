#!/usr/bin/bash

function setup()
{
  module load python/3.8
  source ~/rds/public_databases/software/py38/bin/activate
  for f in .gitignore README.md caprion.ini caprion.ipynb caprion.R docs.sh mkdocs.yml utils
  do cp -r -p ../pilot/${f} pilot; done
}
cp -p ../pilot/autoencoder.md pilot/autoencoder
cp -p ../pilot/gwas2.md pilot/gwas2

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
git add misc/ progs/ peptide_progs/
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
