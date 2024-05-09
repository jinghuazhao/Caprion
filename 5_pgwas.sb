#!/usr/bin/bash

#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --mem=50000
#SBATCH --time=12:00:00

##SBATCH --job-name=_freq
##SBATCH --array=1-22
##SBATCH --mem=50000
##SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/bgen/_freq-0.01_%A_%a.o
##SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/bgen/_freq-0.01_%A_%a.e

##SBATCH --job-name=_fastGWA
##SBATCH --array=1-987
##SBATCH --mem=28800
##SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/pgwas_dr/slurm/_pgwas_fastLR_%A_%a.o
##SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/pgwas_dr/slurm/_pgwas_fastLR_%A_%a.e

#SBATCH --job-name=_qqmanhattan
#SBATCH --array=1-987
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL/qqmanhattanlz/slurm/_qqmanhattan_%A_%a.o
#SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL/qqmanhattanlz/slurm/_qqmanhattan_%A_%a.e

##SBATCH --job-name=_lz
##SBATCH --array=1-987
##SBATCH --mem 6800
##SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL_dr/qqmanhattanlz/slurm/_lz_%A_%a.o
##SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL_dr/qqmanhattanlz/slurm/_lz_%A_%a.e

##SBATCH --job-name=pairs
##SBATCH --array=1-987
##SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL/pairs/slurm/_pairs_%A_%a.o
##SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL/pairs/slurm/_pairs_%A_%a.e

##SBATCH --job-name=miamiplot
##SBATCH --array=1-987
##SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL/miamiplot/slurm/_miamiplot_%A_%a.o
##SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/Caprion_proteomics/analysis/METAL/miamiplot/slurm/_miamiplot_%A_%a.e

export TMPDIR=${HPC_WORK}/work
export caprion=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=
export phenoname=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${analysis}/output/caprion${suffix}.varlist)
export interval=${HPC_WORK}/data/interval

. /etc/profile.d/modules.sh
module purge
export PERL5LIB=
module load rhel8/default-icl
module load samtools/1.13/gcc/zwxn7ug3

function rsid_snpid()
{
  seq 22 | \
  parallel -j1 -C' ' '
    echo {}
    export chr=chr${}
    export snpstats=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/impute_{}_interval.snpstats
    sed '1d' ${snpstats} | cut -f2-6 | awk "
    (\$4<\$5) {
      \$2=\$2+0; print \$2, \$3, \$1, \$2\":\"\$3\"_\"\$4\"_\"\$5
    }{
      \$2=\$2+0; print \$2, \$3, \$1, \$2\":\"\$3\"_\"\$5\"_\"\$4
    } " > ${analysis}/bgen/chr{}.rsid
  '
}

function bgen()
{
  export chr=chr${SLURM_ARRAY_TASK_ID}
  plink2 --bgen ${interval}/${chr}.bgen ref-unknown --sample ${interval}/interval.sample --keep ${analysis}/work/caprion${suffix}.id2 \
         --export bgen-1.2 bits=8 \
         --set-all-var-ids @:#_\$1_\$2 --new-id-max-allele-len 680 \
         --out ${analysis}/bgen/${chr}
  bgenix -g ${analysis}/bgen/${chr}.bgen -index -clobber
}

function freq()
{
  export chr=chr${SLURM_ARRAY_TASK_ID}
  plink2 --bgen ${analysis}/bgen/${chr}.bgen ref-unknown --sample ${analysis}/bgen/${chr}.sample \
         --freq \
         --out ${analysis}/bgen/${chr}-freq
  awk 'NR>1 && $5>=0.001 && $5<=0.999 {print $2}' ${analysis}/bgen/${chr}-freq.afreq > ${analysis}/bgen/${chr}.snplist
}

function freq.01()
{
  export chr=chr${SLURM_ARRAY_TASK_ID}
  plink2 --bgen ${analysis}/bgen/${chr}.bgen ref-unknown --sample ${analysis}/bgen/${chr}.sample \
         --freq \
         --out ${analysis}/bgen/${chr}-freq
  awk 'NR>1 && $5>=0.01 && $5<=0.99 {print $2}' ${analysis}/bgen/${chr}-freq.afreq > ${analysis}/bgen/${chr}-0.01.snplist
}

function fastLR()
# fastGWA without --grm-sparse
{
  export phenocol=${SLURM_ARRAY_TASK_ID}
  export phenoname=$(awk 'NR==ENVIRON["phenocol"]{print $1}' ${analysis}/output/caprion${suffix}.varlist)
  export batch=${1}
  gcta-1.9 --mbgen ${analysis}/bgen/caprion.bgenlist \
           --sample ${analysis}/bgen/caprion.sample \
           --extract ${analysis}/bgen/caprion.snplist \
           --keep ${analysis}/output/caprion${suffix}-${batch}.id \
           --fastGWA-lr \
           --pheno ${analysis}/output/caprion${suffix}-${batch}.pheno --mpheno ${phenocol} \
           --threads 10 \
           --out ${analysis}/pgwas${suffix}/caprion${suffix}-${batch}-${phenoname}

  gcta-1.9 --mbgen ${analysis}/bgen/caprion.bgenlist \
           --sample ${analysis}/bgen/caprion.sample \
           --extract ${analysis}/bgen/caprion.snplist \
           --keep ${analysis}/output/chrX${suffix}-${batch}.id \
           --fastGWA-lr --model-only \
           --pheno ${analysis}/output/caprion${suffix}-${batch}.pheno --mpheno ${phenocol} \
           --threads 10 \
           --out ${analysis}/pgwas${suffix}/caprion${suffix}-${batch}-${phenoname}

  gcta-1.9 --bgen ${analysis}/bgen/chrX.bgen \
           --sample ${analysis}/bgen/chrX.sample \
           --extract ${analysis}/bgen/caprion.snplist --geno 0.1 \
           --keep ${analysis}/output/chrX${suffix}-${batch}.id \
           --load-model ${analysis}/pgwas${suffix}/caprion${suffix}-${batch}-${phenoname}.fastGWA \
           --threads 10 \
           --out ${analysis}/pgwas${suffix}/caprion${suffix}-${batch}-${phenoname}-chrX
  bgzip -f ${analysis}/pgwas${suffix}/caprion${suffix}-${batch}-${phenoname}.fastGWA
  bgzip -f ${analysis}/pgwas${suffix}/caprion${suffix}-${batch}-${phenoname}-chrX.fastGWA
}

# export phenocol=$(grep -n -f ${caprion}/work/caprion.lrlist ${caprion}/work/caprion.varlist | \
#                   tr ':' '\t' | cut -f1 | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]')

function qqmanhattan()
{
  export phenoname=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${analysis}/work/caprion${suffix}.varlist)${suffix}
  gunzip -c ${analysis}/METAL${suffix}/${phenoname}-1.tbl.gz | \
  awk '{if (NR==1) print "chromsome","position","log_pvalue","beta","se";
        else if ($1!=23) print $1,$2,-$12,$10,$11}' | \
  gzip -f > ${analysis}/work/${phenoname}.txt.gz
  R --slave --vanilla --args \
      input_data_path=${analysis}/work/${phenoname}.txt.gz \
      output_data_rootname=${analysis}/METAL${suffix}/qqmanhattanlz/${phenoname}_qq \
      plot_title="${phenoname}" < ~/cambridge-ceu/turboqq/turboqq.r
  if [ ! -f ${analysis}/METAL${suffix}/sentinels/${phenoname}.signals ]; then
     R --slave --vanilla --args \
       input_data_path=${analysis}/work/${phenoname}.txt.gz \
       output_data_rootname=${analysis}/METAL${suffix}/qqmanhattanlz/${phenoname}_manhattan \
       reference_file_path=~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda \
       pvalue_sign=5e-8 \
       plot_title="${phenoname}" < ~/cambridge-ceu/turboman/test.r
  else
    R --slave --vanilla --args \
      input_data_path=${analysis}/work/${phenoname}.txt.gz \
      output_data_rootname=${analysis}/METAL${suffix}/qqmanhattanlz/${phenoname}_manhattan \
      custom_peak_annotation_file_path=${analysis}/METAL${suffix}/vep/${phenoname}.txt \
      reference_file_path=~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda \
      pvalue_sign=5e-8 \
      plot_title="${phenoname}" < ~/cambridge-ceu/turboman/test.r
  fi
  rm ${analysis}/work/${phenoname}.txt.gz
# R --no-save < utils/qqman.R
#   cat <(echo chromosome position) \
#       <(awk 'NR>1{print $1,$2}' ${analysis}/METAL${suffix}/sentinels/${phenoname}.signals) \
#       > ${analysis}/work/${phenoname}.annotate
}

function lz()
{
  module load python/2.7
  export phenoname=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${analysis}/work/caprion${suffix}.varlist)${suffix}
  if [ -f ${analysis}/METAL${suffix}/sentinels/${phenoname}.signals ]; then
     (
       awk 'NR>1{print $6}' ${analysis}/METAL${suffix}/sentinels/${phenoname}.signals | \
       parallel -j1 -C' ' --env analysis --env phenoname '
         zgrep -w {} ${analysis}/METAL${suffix}/${phenoname}-1.tbl.gz | \
         awk -v rsid={} "{print \$1,\$2-5e5,\$2+5e5,rsid}"
       '
     ) | \
     parallel -j1 -C ' ' --env analysis --env phenoname '
     export type=$(awk "\$2==rsid {\$3=\$3 suffix; if(\$3==prot) print \$10}" FS="," rsid={4} prot=${phenoname} suffix=${suffix} \
                   ${analysis}/work/caprion${suffix}.cis.vs.trans)
     if [ {1} != "X" ]; then
       (
         echo -e "Chromosome\tPosition\tMarkerName\tlog10P"
         gunzip -c ${analysis}/METAL${suffix}/${phenoname}-1.tbl.gz | \
         awk -v chr={1} -v start={2} -v end={3} -v OFS="\t" "\$1==chr && \$2>=start && \$2<end {print \$1,\$2,\$3,-\$12}" | \
         sort -k1,1n -k2,2n
       ) > ${analysis}/work/${phenoname}-{4}.lz
       locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${analysis}/work/${phenoname}-{4}.lz \
                 --delim tab title="${phenoname}-{4} ($type)" \
                 --markercol MarkerName --pvalcol log10P --no-transform --chr {1} --start {2} --end {3} --cache None \
                 --no-date --plotonly --prefix=${phenoname} --rundir ${analysis}/METAL${suffix}/qqmanhattanlz --refsnp {4}
     else
       (
         echo -e "Chromosome\tPosition\tMarkerName\tlog10P"
         gunzip -c ${analysis}/METAL${suffix}/${phenoname}-1.tbl.gz | \
         awk -v chr={1} -v start={2} -v end={3} -v OFS="\t" "\$1==chr && \$2>=start && \$2<end {print \$1,\$2,\$3,-\$12}}" | \
         sort -k1,1n -k2,2n | \
         sed "s/_[A-Z]*_[A-Z]*//" | cut -f1-3,12 | sed "s/X/chr23/"
       ) > ${analysis}/work/${phenoname}-chrX-{4}.lz
       locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${analysis}/work/${phenoname}-chrX-{4}.lz \
                 --delim tab title="${phenoname}-chr{4} ($type)" \
                 --markercol MarkerName --pvalcol log10P --no-transform --chr {1} --start {2} --end {3} --cache None \
                 --no-date --plotonly --prefix=${phenoname}-chrX --rundir ${analysis}/METAL${suffix}/qqmanhattanlz \
                 --refsnp $(echo chr{4} | sed "s/_[A-Z]*_[A-Z]*//")
       rm ${analysis}/work/${phenoname}-{4}.lz
     fi
     '
  fi
}
# HSPB1_dr_rs114800762.pdf would fail with 1000G_Nov2014 --hg19 --pop EUR, so we get around (actually copy) with interval/interval37
# #SBATCH --array=462
# Warning: rs114800762 is not the current name in genome build (should be: rs3982708)
# export old=rs114800762
# export new=rs3982708
# On the other hand, the following list is missing from interval/interval37 which are copied from 1000G_Nov2014/hg19
# cp lz-interval37/*pdf .
# diff <(ls lz-interval37) <(ls lz) | grep '>' | sed 's/> //' | xargs -I {} cp lz/{} .
# ACE_dr_rs587729126.pdf
# B3GN7_dr_rs537179722.pdf
# CFAB_dr_rs575150094.pdf
# CYTC_dr_rs55727201.pdf
# FCN3_dr_rs75318474.pdf
# I13R1_dr_rs55727201.pdf
# ICAM1_dr_rs587745765.pdf
# ICAM2_dr_rs587745765.pdf
# LFA3_dr_rs587729126.pdf
# NQO2_dr_rs370546929.pdf
# PODXL_dr_rs587729126.pdf

#echo "AMY1A; AMY1B; AMY1C" | tr ';' '\n' | sed 's/ //' | grep -f - ~/INF/csd3/glist-hg19
#1 104198140 104207173 AMY1A
#1 104230039 104239073 AMY1A
#1 104292278 104301311 AMY1A
#1 104198324 104207172 AMY1B
#1 104230040 104238889 AMY1B
#1 104292462 104301310 AMY1B
#1 104198302 104207172 AMY1C
#1 104230040 104238911 AMY1C
#1 104292440 104301310 AMY1C

#echo "APOC2" | tr ';' '\n' | sed 's/ //' | grep -f - ~/INF/csd3/glist-hg19
#19 45449238 45452822 APOC2
#19 45445494 45452822 APOC4-APOC2

#echo "APOC4" | tr ';' '\n' | sed 's/ //' | grep -f - ~/INF/csd3/glist-hg19
#19 45445494 45448753 APOC4
#19 45445494 45452822 APOC4-APOC2

#echo "C4B; C4B_2" | tr ';' '\n' | sed 's/ //' | grep -f - -w ~/INF/csd3/glist-hg19
#6 31949833 31970458 C4B
#6 31982571 32003195 C4B
#6 31949833 31970458 C4B_2
#6 31982571 32003195 C4B_2

#echo "HIST1H4A; HIST1H4B; HIST1H4C; HIST1H4D; HIST1H4E; HIST1H4F; HIST1H4H; HIST1H4I; HIST1H4J; HIST1H4K; HIST1H4L; HIST2H4A; HIST2H4B; HIST4H4" | tr ';' '\n' | sed 's/ //' | grep -f - ~/INF/csd3/glist-hg19
#6 26021906 26022278 HIST1H4A
#6 26027123 26027480 HIST1H4B
#6 26104175 26104565 HIST1H4C
#6 26188937 26189304 HIST1H4D
#6 26204872 26205249 HIST1H4E
#6 26240653 26241021 HIST1H4F
#6 26285353 26285727 HIST1H4H
#6 27107087 27107457 HIST1H4I
#6 27791902 27792258 HIST1H4J
#6 27798951 27799305 HIST1H4K
#6 27840925 27841289 HIST1H4L
#1 149804220 149804616 HIST2H4A
#1 149832329 149832725 HIST2H4A
#1 149804220 149804616 HIST2H4B
#1 149832329 149832725 HIST2H4B
#12 14923653 14924065 HIST4H4

#echo "HBA1; HBA2" | tr ';' '\n' | sed 's/ //' | grep -f - ~/INF/csd3/glist-hg19
#16 226678 227520 HBA1
#16 222845 223709 HBA2

function pairs()
# all variants
{
  export caprion=~/Caprion
  export analysis=${caprion}/analysis
  export pgwas=${analysis}/pgwas
  export pilot=${caprion}/pilot
  export work=${analysis}/work
  export TMPDIR=${HPC_WORK}/work
  export uniprot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' ~/Caprion/analysis/work/caprion.list)
  export protx=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ~/Caprion/analysis/work/caprion.list | sed 's/^\([0-9].*\)/X\1/g')
  export prot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ~/Caprion/analysis/work/caprion.list)

  gunzip -c ${pilot}/bgen/${uniprot}_invn-plink2.gz | cut -f1-3,6,9,10,12 | gzip -f > ${work}/${prot}-pilot-1.gz
  gunzip -c ${pilot}/bgen2/${protx}_All_invn-plink2.gz | cut -f1-3,6,9,10,12 | gzip -f > ${work}/${prot}-pilot-2.gz
  gunzip -c ${pilot}/bgen3/protein-${prot}_invn-plink2.gz | cut -f1-3,6,9,10,12 | gzip -f > ${work}/${prot}-pilot-3.gz
  seq 3 | parallel -C' ' 'gunzip -c ${pgwas}/caprion-{}-${prot}.fastGWA.gz | cut -f1-3,4,8-10 | gzip -f> ${work}/${prot}-{}.gz'

  Rscript -e '
    suppressMessages(library(dplyr))
    prot <- Sys.getenv("prot")
    work <- Sys.getenv("work")
    b1 <- read.delim(file.path(work,paste0(prot,"-pilot-1.gz"))) %>% rename(CHR=X.CHROM,SNP=ID,b1=BETA,a1=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b2 <- read.delim(file.path(work,paste0(prot,"-pilot-2.gz"))) %>% rename(CHR=X.CHROM,SNP=ID,b2=BETA,a2=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b3 <- read.delim(file.path(work,paste0(prot,"-pilot-3.gz"))) %>% rename(CHR=X.CHROM,SNP=ID,b3=BETA,a3=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b <- b1 %>% left_join(b2) %>% left_join(b3) %>% mutate(b2=if_else(a1==a2,b2,-b2),b3=if_else(a1==a3,b3,-b3))
    cor(b[c("b1","b2","b3")],use="complete.obs",method="pearson")
    png(file.path("work",paste0(prot,"-pilot.png")),width=10,height=8,units="in",res=300)
    pairs(b[paste0("b",1:3)],cex=0.5,pch=19,upper.panel=NULL)
    dev.off()
    b1 <- read.delim(file.path(work,paste0(prot,"-1.gz"))) %>% rename(b1=BETA,a1=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b2 <- read.delim(file.path(work,paste0(prot,"-2.gz"))) %>% rename(b2=BETA,a2=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b3 <- read.delim(file.path(work,paste0(prot,"-3.gz"))) %>% rename(b3=BETA,a3=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b <- b1 %>% left_join(b2) %>% left_join(b3) %>% mutate(b2=if_else(a1==a2,b2,-b2),b3=if_else(a1==a3,b3,-b3))
    cor(b[c("b1","b2","b3")],use="complete.obs",method="pearson")
    png(file.path("work",paste0(prot,".png")),width=10,height=8,units="in",res=300)
    pairs(b[paste0("b",1:3)],cex=0.5,pch=19,upper.panel=NULL)
    dev.off()
    p1 <- read.delim(file.path(work,paste0(prot,"-pilot-1.gz"))) %>% rename(CHR=X.CHROM,SNP=ID,p1=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p1=-log10(p1))
    p2 <- read.delim(file.path(work,paste0(prot,"-pilot-2.gz"))) %>% rename(CHR=X.CHROM,SNP=ID,p2=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p2=-log10(p2))
    p3 <- read.delim(file.path(work,paste0(prot,"-pilot-3.gz"))) %>% rename(CHR=X.CHROM,SNP=ID,p3=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p3=-log10(p3))
    p <- p1 %>% left_join(p2) %>% left_join(p3)
    cor(p[c("p1","p2","p3")],use="complete.obs",method="pearson")
    png(file.path("work",paste0(prot,"-pilot-p.png")),width=10,height=8,units="in",res=300)
    pairs(p[paste0("p",1:3)],cex=0.5,pch=19,upper.panel=NULL)
    dev.off()
    p1 <- read.delim(file.path(work,paste0(prot,"-1.gz"))) %>% rename(p1=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p1=-log10(p1))
    p2 <- read.delim(file.path(work,paste0(prot,"-2.gz"))) %>% rename(p2=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p2=-log10(p2))
    p3 <- read.delim(file.path(work,paste0(prot,"-3.gz"))) %>% rename(p3=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p3=-log10(p3))
    p <- p1 %>% left_join(p2) %>% left_join(p3)
    cor(p[c("p1","p2","p3")],use="complete.obs",method="pearson")
    png(file.path("work",paste0(prot,"-p.png")),width=10,height=8,units="in",res=300)
    pairs(p[paste0("p",1:3)],cex=0.5,pch=19,upper.panel=NULL)
    dev.off()
  '
  rm ${work}/${prot}-pilot-?.gz ${work}/${prot}-?.gz
}

function pairs_p()
# variants below 1e-5 cutoff
{
  export caprion=~/Caprion
  export analysis=${caprion}/analysis
  export pgwas=${analysis}/pgwas
  export pilot=${caprion}/pilot
  export work=${analysis}/work
  export TMPDIR=${HPC_WORK}/work
  export uniprot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' ~/Caprion/analysis/work/caprion.list)
  export protx=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ~/Caprion/analysis/work/caprion.list | sed 's/^\([0-9].*\)/X\1/g')
  export prot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ~/Caprion/analysis/work/caprion.list)

  seq 3 | parallel -C' ' 'gunzip -c ${pgwas}/caprion-{}-${prot}.fastGWA.gz | awk "NR==1||\$10<=1e-5" | cut -f1-3,4,8-10 | gzip -f > ${work}/${prot}-{}.gz'

  Rscript -e '
    suppressMessages(library(dplyr))
    prot <- Sys.getenv("prot")
    work <- Sys.getenv("work")
    b1 <- read.delim(file.path(work,paste0(prot,"-1.gz"))) %>% rename(b1=BETA,a1=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b2 <- read.delim(file.path(work,paste0(prot,"-2.gz"))) %>% rename(b2=BETA,a2=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b3 <- read.delim(file.path(work,paste0(prot,"-3.gz"))) %>% rename(b3=BETA,a3=A1) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-SE,-P)
    b <- b1 %>% left_join(b2) %>% left_join(b3) %>% mutate(b2=if_else(a1==a2,b2,-b2),b3=if_else(a1==a3,b3,-b3))
    cor(b[c("b1","b2","b3")],use="complete.obs",method="pearson")
    png(file.path("~/Caprion/analysis/METAL/pairs",paste0(prot,".png")),width=10,height=8,units="in",res=300)
    pairs(b[paste0("b",1:3)],cex=0.5,pch=19,upper.panel=NULL)
    dev.off()
    p1 <- read.delim(file.path(work,paste0(prot,"-1.gz"))) %>% rename(p1=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p1=-log10(p1))
    p2 <- read.delim(file.path(work,paste0(prot,"-2.gz"))) %>% rename(p2=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p2=-log10(p2))
    p3 <- read.delim(file.path(work,paste0(prot,"-3.gz"))) %>% rename(p3=P) %>%
          mutate(chrpos=paste0(CHR,":",POS)) %>% select(-CHR,-POS,-SNP,-BETA,-SE) %>% mutate(p3=-log10(p3))
    p <- p1 %>% left_join(p2) %>% left_join(p3)
    cor(p[c("p1","p2","p3")],use="complete.obs",method="pearson")
    png(file.path("~/Caprion/analysis/METAL/pairs",paste0(prot,"-p.png")),width=10,height=8,units="in",res=300)
    pairs(p[paste0("p",1:3)],cex=0.5,pch=19,upper.panel=NULL)
    dev.off()
  '
  rm ${work}/${prot}-?.gz
}

function miamiplot()
{
  export TMPDIR=/rds/user/jhz22/hpc-work/work
  export caprion=~/rds/projects/Caprion_proteomics/pilot
  export prot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${caprion}/2020.id | sed -r 's/^X([0-9])/\1/')
  export uniprot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' ${caprion}/2020.id)
  (
    echo snpid chr pos rsid p z pr zr
    join \
    <(gunzip -c ~/Caprion/analysis/pgwas/caprion-1-${prot}.fastGWA.gz | \
      awk 'NR>1{if($4<$5) {a1=$4;a2=$5} else {a1=$5;a2=$4}; snpid=$1":"$2"_"a1"_"a2; print snpid,$1,$3,$2,$10,$8/$9}' | \
      sort -k1,1 \
     ) \
    <(gunzip -c ~/Caprion/analysis/pgwas/caprion-2-${prot}.fastGWA.gz  | \
      awk 'NR>1{if($4<$5) {a1=$4;a2=$5} else {a1=$5;a2=$4}; snpid=$1":"$2"_"a1"_"a2; print snpid,$10,$8/$9}' | \
      sort -k1,1 \
     )
  ) | gzip -f > ${analysis}/work/${uniprot}-${prot}-phase1-phase2.dat.gz
  module load gcc/6
  Rscript -e '
    analysis <- Sys.getenv("analysis")
    prot <- Sys.getenv("prot")
    uniprot <- Sys.getenv("uniprot")
    uniprot_prot <- read.table(file.path(analysis,"work",paste(uniprot,prot,"phase1-phase2.dat.gz",sep="-")),as.is=TRUE,header=TRUE)
    png(file.path(analysis,"METAL","miamiplot",paste(uniprot,prot,"phase1-phase2-fastGWA.png",sep="-")),res=300,width=12,height=10,units="in")
    gap::miamiplot(uniprot_prot,chr="chr",bp="pos",p="pr",pr="p",snp="rsid",cex=0.4,ylab="Phase II (top)/I (bottom)")
    z <- function(p) qnorm(p/2,lower.tail=FALSE)
    cor_z1_z2 <- with(uniprot_prot,cor(z,zr,use="complete.obs",method="pearson"))
    cor_p1_p2 <- with(uniprot_prot,cor(sign(z)*z(p),sign(zr)*z(pr),use="complete.obs",method="pearson"))
    sign_test_p1_p2 <- with(uniprot_prot,wilcox.test(p, pr, alternative="less", paired=TRUE, na.action="na.omit"))
    legend("topright",sprintf("Pearson r(z1,z2)=%.4f, r(z(p1),z(p2))=%.4f,Wilcoxon signed rank test p=%.4f",
           cor_z1_z2,cor_p1_p2,with(sign_test_p1_p2,p.value)))
    dev.off()
  '
  rm ${analysis}/work/${uniprot}-${prot}-phase1-phase2.dat.gz
}

function miamiplot2()
{
  export TMPDIR=/rds/user/jhz22/hpc-work/work
  export caprion=~/rds/projects/Caprion_proteomics/pilot
  export analysis=~/Caprion/analysis
  export prot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${caprion}/2020.id | sed -r 's/^X([0-9])/\1/')
  export uniprot=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' ${caprion}/2020.id)
  module load gcc/6
  if [ ! -f ${analysis}/work/${uniprot}-${prot}-phase1.dat.gz ]; then
     (
       echo snpid chr pos rsid p z
       gunzip -c ~/Caprion/analysis/pgwas/caprion-1-${prot}.fastGWA.gz | \
       awk 'NR>1{if($4<$5) {a1=$4;a2=$5} else {a1=$5;a2=$4}; snpid=$1":"$3"_"a1"_"a2; print snpid,$1,$3,$2,$10,$8/$9}' | \
       sort -k2,2n -k3,3n
     ) | \
     gzip -f > ${analysis}/work/${uniprot}-${prot}-phase1.dat.gz
  fi
  if [ ! -f ${analysis}/work/${uniprot}-${prot}-phase2.dat.gz ]; then
     (
       echo snpid chr pos rsid p z
       gunzip -c ~/Caprion/analysis/pgwas/caprion-2-${prot}.fastGWA.gz  | \
       awk 'NR>1{if($4<$5) {a1=$4;a2=$5} else {a1=$5;a2=$4}; snpid=$1":"$3"_"a1"_"a2; print snpid,$1,$3,$2,$10,$8/$9}' | \
       sort -k2,2n -k3,3n
     ) | \
     gzip -f > ${analysis}/work/${uniprot}-${prot}-phase2.dat.gz
  fi
  Rscript -e '
  caprion <- Sys.getenv("caprion");analysis <- Sys.getenv("analysis"); prot <- Sys.getenv("prot"); uniprot <- Sys.getenv("uniprot")
  protein <- prot
  suppressMessages(require(dplyr))
  suppressMessages(require(gap))
  cvt <- read.csv("~/Caprion/analysis/work/caprion.cis.vs.trans")
  annot <- subset(cvt,prot==protein) %>%
           mutate(name=if_else(cis,paste(Gene,SNP,sep="-"),SNP))
  gwas1 <- read.table(file.path(analysis,"work",paste(uniprot,prot,"phase2.dat.gz",sep="-")),as.is=TRUE,header=TRUE) %>% filter(!is.na(p))
  gwas2 <- read.table(file.path(analysis,"work",paste(uniprot,prot,"phase1.dat.gz",sep="-")),as.is=TRUE,header=TRUE) %>% filter(!is.na(p))
  png(file.path(analysis,"METAL","miamiplot",paste(uniprot,prot,"phase1-phase2-fastGWA.png",sep="-")),res=300,width=12,height=10,units="in")
  chrmaxpos <- miamiplot2(gwas1,gwas2,name1="Batch 2",name2="Batch 1",z1="z",z2="z") %>%
               filter(chr!=-Inf & maxpos!=-Inf & genomestartpos!=-Inf & labpos!=-Inf)
  labelManhattan(chr=annot$SNPChrom,pos=annot$SNPPos,name=annot$name,gwas1,gwasZLab="z",chrmaxpos=chrmaxpos)
  dev.off()
  '
  rm ${analysis}/work/${uniprot}-${prot}-phase?.dat.gz
}

#bgen
#freq
#freq.01
#fastLR 1
#fastLR 2
#fastLR 3
qqmanhattan
#lz
#pairs_p
#miamiplot
#miamiplot2

# --- legacy code ---

function fastGWA()
# fastGWA mixed model
{
  gcta-1.9 --mbgen ${caprion}/bgen/caprion.bgenlist --grm-sparse ${caprion}/output/caprion-spgrm \
           --sample ${caprion}/bgen/caprion.sample \
           --fastGWA-mlm --pheno ${caprion}/output/caprion.pheno --mpheno ${SLURM_ARRAY_TASK_ID} --threads 10 \
           --out ${caprion}/work/caprion-${phenoname}

  gcta-1.9 --mbgen ${caprion}/bgen/caprion.bgenlist --grm-sparse ${caprion}/output/caprion-spgrm \
           --sample ${caprion}/bgen/caprion.sample \
           --fastGWA-mlm --model-only --pheno ${caprion}/output/caprion.pheno --mpheno ${SLURM_ARRAY_TASK_ID} \
           --keep ${caprion}/output/chrX.idlist --threads 10 \
           --out ${caprion}/work/caprion-${phenoname}

  gcta-1.9 --bgen ${caprion}/bgen/chrX.bgen \
           --sample ${caprion}/bgen/caprion.sample \
           --load-model ${caprion}/work/caprion-${phenoname}.fastGWA \
           --extract ${caprion}/output/chrX.snplist --geno 0.1 --threads 10 \
           --out ${caprion}/work/caprion-${phenoname}-chrX
  bgzip -f ${caprion}/work/caprion-${phenoname}*fastGWA
}
