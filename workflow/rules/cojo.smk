'''Snakefile for conditional analysis'''

import os
import pandas
from datetime import datetime
from datetime import date

RWD = os.getcwd()
today_date = date.today()
shell.prefix('module load gcc/6 plink/2.00-alpha python/3.7 ceuadmin/stata texlive')

m = pandas.read_table("work/caprion.merge").drop(['FreqSE','MinFreq','MaxFreq','Direction','HetISq','HetChiSq','HetDf','logHetP'], axis=1)
m = m[m.prot=="ERAP2"]
sentinel = m.assign(id=m.prot+","+m.MarkerName+","+m.Chromosome.astype(str)+","+m.Position.astype(str)).set_index(['prot','MarkerName'], drop=False)

rule all:
    input:
        expand("results/{prot}-{rsid}-{chr}:{pos}",
               zip,
               prot=sentinel['prot'],
               rsid=sentinel['MarkerName'],
               chr=sentinel['Chromosome'],
               pos=sentinel['Position'])

rule cojo:
    params:
        bfile="~/Caprion/pilot/data/merged_imputation",
        prot="{prot}",
        rsid="{rsid}",
        chr="{chr}",
        pos="{pos}"
    output:
        "results/{prot}-{rsid}-{chr}:{pos}"
    shell:
        "../scripts/cojo.sh {params.bfile} {params.prot} {params.rsid} {params.chr} {params.pos}"

# https://tinyheero.github.io/2019/08/30/wildcards-in-snakemake.html
