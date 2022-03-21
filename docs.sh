#!/usr/bin/bash

function setup()
{
  module load python/3.7
  source ~/COVID-19/py37/bin/activate
}

setup
mkdocs build
mkdocs gh-deploy

git add .gitignore
git commit -m ".gitignore"
git add README.md
git commit -m "README"
git add 1_pca_projection.sh 2_ggm.R 3_wgcna.R 4_pca_clustering.R
git commit -m "Analysis"
git add caprion.ini caprion.R caprion.ipynb
git commit -m "Primary Code"
git add gwas2.md
git commit -m "gwas2"
git add tests.md
git commit -m "tests"
git add utils docs.sh
git commit -m "utils"
git add mkdocs.yml
git commit -m "mkdocs.yml"
git push
# git remote set-url origin https://github.com/jinghuazhao/Caprion.git
