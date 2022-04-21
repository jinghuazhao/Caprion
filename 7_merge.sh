#!/usr/bin/bash

export caprion=~/Caprion/analysis
export TMPDIR=${HPC_WORK}/work

function setup()
{
  if [ -d ${caprion}/METAL/sentinels ]; then mkdir -p ${caprion}/METAL/sentinels
}

setup
