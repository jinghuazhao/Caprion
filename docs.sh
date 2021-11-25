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
git add caprion.ini caprion.R caprion.ipynb
git commit -m "Primary Code"
git add gwas2.md
git commit -m "gwas2"
git add tests.md
git commit -m "tests"
git add mkdocs.yml
git commit -m "mkdocs.yml"
git push
# git remote set-url origin https://github.com/jinghuazhao/Caprion.git
