#!/usr/bin/bash

git add README.md
git commit -m "README"
git add caprion.ini caprion.R caprion.ipynb
git commit -m "Primary Code"
git add utils
git commit -m "utilities"
git push
