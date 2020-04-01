#!/usr/bin/bash

export rt="caprion-invn"
cut -d' ' -f6 ${rt}.sentinels | sed '1d' | sort | uniq > ${rt}.snp

