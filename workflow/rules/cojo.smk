import pandas
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
        bfile="data/merged_imputation",
        prot="{prot}",
        rsid="{rsid}",
        chr="{chr}",
        pos="{pos}",
        chunksize=10
    output:
        "results/{prot}-{rsid}-{chr}:{pos}"
    shell:
        "workflow/scripts/cojo.sh {params.bfile} {params.prot} {params.rsid} {params.chr} {params.pos}"
