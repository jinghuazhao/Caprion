#!/usr/bin/bash

function setup()
{
module load python/3.7
source ~/COVID-19/py37/bin/activate
}

# setup
mkdocs build
mkdocs gh-deploy

git remote set-url origin https://github.com/jinghuazhao/Caprion.git
git add .gitignore
git commit -m ".gitignore"
git add README.md
git commit -m "README"
git add caprion.ini caprion.R caprion.ipynb
git commit -m "Primary Code"
git add utils
git commit -m "utilities"
git add mkdocs.yml
git commit -m "mkdocs.ytml"
git push
