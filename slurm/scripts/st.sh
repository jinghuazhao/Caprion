# --unlock to be used optionally
snakemake -c --profile cojo -s cojo/rules/cojo.yaml
snakemake -c --cluster-config cojo/config.yaml -s cojo/rules/cojo.yaml -d cojo
