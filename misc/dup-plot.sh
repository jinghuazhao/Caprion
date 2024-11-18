#!/usr/bin/bash

function step3_pqtl_summary()
{
cat <<'EOL'> ${root}/dup-step3.sh
#!/usr/bin/bash

#SBATCH --job-name=_3-LABEL
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem

#SBATCH --export ALL
#SBATCH --array=_array_
#SBATCH --output=ROOT/sentinels/slurm/_step3_%A_%a.o
#SBATCH --error=ROOT/sentinels/slurm/_step3_%A_%a.e

export TMPDIR=${HPC_WORK}/work
export isotope=$(head -1 ${root}/ZWK.pheno | awk -v n=${SLURM_ARRAY_TASK_ID} '{print $(n+2)}')
export PERL5LIB=

export root=${root}
export dir=${root}/qqmanhattanlz
export isotope=${1}

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/R
module load samtools/1.13/gcc/zwxn7ug3
module load perl/5.26.3_system/gcc-8.4.1-4cl2czq
module load libiconv/1.16/intel/64iicvbf

function qqmanhattan()
{
  module load python/3.7
  source ~/COVID-19/py37/bin/activate
# head -1 ${root}/ZWK.pheno | cut -d' ' -f1-2 --complement | cut -d' ' -f${SLURM_ARRAY_TASK_ID} | tr ' ' '\n' | \
  echo ${isotope} | \
  parallel -C' ' --env root '
  (
    echo chromosome position pvalue beta se
    gunzip -c ${root}/ZWK-1-{}.fastGWA.gz | \
    awk "NR>1{print \$1,\$3,\$10,\$8,\$9}" | \
    sort -k1,1n -k2,2n
  ) > ${root}/work/{}.txt
  R --slave --vanilla --args \
      input_data_path=${root}/work/{}.txt \
      output_data_rootname=${dir}/{}_qq \
      plot_title="{}" < ~/cambridge-ceu/turboqq/turboqq.r
  if [ ! -f ${root}/sentinels/{}.signals ]; then
     R --slave --vanilla --args \
       input_data_path=${root}/work/{}.txt \
       output_data_rootname=${dir}/{}_manhattan \
       reference_file_path=~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda \
       pvalue_sign=5e-8 \
       plot_title="{}" < ~/cambridge-ceu/turboman/test.r
  fi
  rm ${root}/work/{}.txt
  if [ -f ${dir}/{}_manhattan.png ]; then
     convert +append ${dir}/{}_manhattan.png ${dir}/{}_qq.png -resize x500 -density 300 ${dir}/{}_qqmanhattan.png
     convert ${dir}/{}_qqmanhattan.png -quality 0 ${dir}/{}_qqmanhattan.jp2
     img2pdf -o ${dir}/{}_qqmanhattan.pdf ${dir}/{}_qqmanhattan.jp2
     rm ${dir}/{}_qqmanhattan.jp2
  fi
  '
  deactivate
}

function lz_autosomes()
{
  module load python/2.7
  (
#   awk '$1==ENVIRON["isotope"] && $2!=23 {print $6, $7}' ${root}/${isotope}.signals | \
    awk '$5==ENVIRON["isotope"] && $1!=23 {print $5, $6}' ${root}/dup.signals | \
    parallel -j1 -C' ' --env root '
      zgrep -w {2} ${root}/ZWK-1-{1}.fastGWA.gz | \
      awk -v isotope={1} -v rsid={2} "{print \$1,\$3-5e5,\$3+5e5,isotope,rsid}"
    '
  ) | \
  parallel -j1 -C ' ' --env root '
    (
      gunzip -c ${root}/ZWK-1-{4}.fastGWA.gz | \
      awk -v OFS="\t" "NR==1 {print \$1,\$3,\$10,\$2}"
      gunzip -c ${root}/ZWK-1-{4}.fastGWA.gz | \
      awk -v chr={1} -v start={2} -v end={3} -v OFS="\t" "\$1==chr && \$3>=start && \$3<end {print \$1,\$3,\$10,\$2}" | \
      sort -k1,1n -k2,2n
    ) > ${root}/work/{4}-{5}.lz
    locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${root}/work/{4}-{5}.lz \
              --delim tab title="{4}-{5}" \
              --markercol SNP --pvalcol P --chr {1} --start {2} --end {3} --cache None \
              --no-date --plotonly --prefix={4} --rundir ${root}/qqmanhattanlz --refsnp {5}
    rm ${root}/work/{4}-{5}.lz
  '
  module unload python/2.7
}

function lz_X()
{
  module load python/2.7
  (
  # awk '$1==ENVIRON["isotope"] && $2==23 {print $6, $7}' ${root}/${isotope}.signals | \
    awk '$5==ENVIRON["isotope"] && $1==23 {print $5, $6}' ${root}/dup.signals | \
    parallel -j1 -C' ' --env root '
      zgrep -w {2} ${root}/ZWK-1-{1}-chrX.fastGWA.gz | \
      awk -v isotope={1} -v rsid={2} "{print \$1,\$3-5e5,\$3+5e5,isotope,rsid}" | sed "s/X/chr23/;s/_[A-Z]*_[A-Z]*//"
    '
  ) | \
  parallel -j1 -C ' ' --env root '
    (
      gunzip -c ${root}/ZWK-1-{4}-chrX.fastGWA.gz | \
      awk -v OFS="\t" "NR==1 {print \$1,\$3,\$10,\$2}"
      gunzip -c ${root}/ZWK-1-{4}-chrX.fastGWA.gz | \
      awk -v chr={1} -v start={2} -v end={3} -v OFS="\t" "\$1==chr && \$3>=start && \$3<end {print \$1,\$3,\$10,\$2}" | \
      sort -k1,1n -k2,2n | \
      sed "s/X/chr23/;s/_[A-Z]*_[A-Z]*//"
    ) > ${root}/work/{4}-{5}.lz
    locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${root}/work/{4}-{5}.lz \
              --delim tab title="{4}-{5}" \
              --markercol SNP --pvalcol P --chr {1} --start {2} --end {3} --cache None \
              --no-date --plotonly --prefix={4} --rundir ${root}/qqmanhattanlz \
              --refsnp {5}
    rm ${root}/work/{4}-{5}.lz
  '
  module unload python/2.7
}

function mean_by_genotype()
{
  for batch in 1
  do
    export batch=${batch}
    export out=${root}/means/${batch}-${isotope}-${pqtl}
    if [ ! -f ${out}.dat ]; then
       plink-2 --bgen ${analysis}/bgen/chr${chr}.bgen ref-unknown \
               --sample ${analysis}/bgen/chr${chr}.sample \
               --chr ${chr} --from-bp ${bp} --to-bp ${bp} \
               --keep ${analysis}/output/caprion-${batch}.id \
               --pheno ${root}/ZWK.pheno --pheno-name ${isotope} \
               --recode oxford \
               --out ${out}
       paste <(awk 'NR>2{print $1,$5}' ${out}.sample) \
             <(awk '{for(i=0;i<(NF-5)/3;i++) print $1,$2,$3,$4,$5, $(6+i),$(7+i),$(8+i)}' ${out}.gen) > ${out}.dat
       rm ${out}.gen ${out}.sample ${out}.log
    fi
  done
  Rscript -e '
     options(width=120)
     root <- Sys.getenv("root")
     isotope <- Sys.getenv("isotope")
     pqtl <- Sys.getenv("pqtl")
     invisible(suppressMessages(sapply(c("dplyr","ggplot2","ggpubr"),require,character.only=TRUE)))
     process_batch <- function(batch, genotypes=c("100","010","001"))
     {
       datfile <- file.path(root,"means",paste(batch,isotope,pqtl,sep="-"))
       dat <- read.table(paste0(datfile,".dat"),
                         colClasses=c("character","numeric","character","character","integer","character","character",rep("numeric",3)),
                         col.names=c("IID","Phenotype","chr","rsid","pos","A1","A2","g1","g2","g3")) %>%
              mutate(g=paste0(round(g1),round(g2),round(g3)),
                     Genotype=as.factor(case_when(g == genotypes[1] ~ paste0(A1,"/",A1),
                                                  g == genotypes[2] ~ paste0(A1,"/",A2),
                                                  g == genotypes[3] ~ paste0(A2,"/",A2),
                                                  TRUE ~ "---"))) %>%
              filter(Genotype!="---")
       means <- group_by(dat,Genotype) %>%
                summarise(N=n(),Mean=mean(Phenotype))
       invisible(list(dat=dat,means=means))
     }
     v <- m <- list()
     for (batch in 1)
     {
         x <- process_batch(batch)
         v[[batch]] <- ggplot(with(x,dat), aes(x=Genotype, y=Phenotype, fill=Genotype)) +
                       geom_violin() +
                       geom_boxplot(width=0.1) +
                       theme_minimal()
         m[[batch]] <- ggtexttable(with(x,means), rows = NULL, theme = ttheme("mOrange"))
     }
     p <- ggarrange(v[[1]],m[[1]]],ncol=1,nrow=2,labels=c("1. ZWK"))
     ggsave(file.path(root,"means",paste0(isotope,"-",pqtl,"-genotype.png")),device="png",width=16, height=10, units="in")
  '
}

export -f mean_by_genotype

awk '$1==ENVIRON["isotope"]{gsub(/23/,"X",$1);print $1,$3,$2}' ${root}/dup.merge | \
parallel -C' ' '
  export chr={1}
  export bp={2}
  export pqtl={3}
  mean_by_genotype
'

for cmd in qqmanhattan lz_autosomes lz_X; do $cmd; done
EOL

sed -i "s|ROOT|${root}|;s|_array_|${array}|" ${root}/dup-step3.sh
chmod +x ${root}/dup-step3.sh
}

export analysis=~/Caprion/analysis
export root=~/Caprion/analysis/dup
export pheno=${root}/ZWK.pheno
step3_pqtl_summary

isotopes=$(head -1 ${root}/ZWK.pheno | awk '{for(i=3;i<=NF;i++) print $i}')
echo "$isotopes" | parallel -j 10 --bar ${root}/dup-step3.sh
