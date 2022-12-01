#!/usr/bin/bash

for f in $(ls ); do diff $f ../../pilot/utils/$f; done

rsync -avP . ../../pilot/utils
