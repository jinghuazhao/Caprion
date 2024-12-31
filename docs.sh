#!/usr/bin/bash

function setup()
{
  module load python/3.8
  source ~/rds/public_databases/software/py38/bin/activate
}

if [ "$(uname -n | sed 's/-[0-9]*$//')" == "login-q" ]; then
   module load ceuadmin/libssh/0.10.6-icelake
   module load ceuadmin/openssh/9.7p1-icelake
fi

setup
mkdocs build
mkdocs gh-deploy

git add index.html
git commit -m 'LocusZoom.js plots'
git add workflow
git commit -m "Snakemake profiles"
git add .gitignore
git commit -m ".gitignore"
git add README.md
git commit -m "README"
git add docs.sh
git commit -m "docs.sh"
git add docs/
git commit -m "docs"
git add pilot/
git commit -m "Pilot studies"
git add misc/ progs/ peptide_progs/
git commit -m "Analysis"
git add mkdocs.yml
git commit -m "mkdocs.yml"
git push
# git remote set-url origin https://github.com/jinghuazhao/Caprion.git
