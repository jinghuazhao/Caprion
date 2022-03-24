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
git add README.md analysis/README.md pilot/README.md
git commit -m "README"
git add analysis/1_pca_projection.sh analysis/2_ggm.R analysis/3_wgcna.R analysis/4_pca_clustering.R analysis/5_pgwas.sh analysis/6_meta_analysis.sh
git commit -m "Analysis"
git add pilot/caprion.ini pilot/caprion.R pilot/caprion.ipynb
git commit -m "Primary Code"
git add pilot/gwas2.md
git commit -m "gwas2"
git add pilot/tests.md
git commit -m "tests"
git add pilot/utils pilot/docs.sh
git commit -m "utils"
git add mkdocs.yml pilot/mkdocs.yml
git commit -m "mkdocs.yml"
git push
# git remote set-url origin https://github.com/jinghuazhao/Caprion.git
