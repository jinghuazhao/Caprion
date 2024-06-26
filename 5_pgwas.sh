#!/usr/bin/bash

export caprion=~/Caprion/pilot
export analysis=~/Caprion/analysis
export suffix=_dr
export interval=${HPC_WORK}/data/interval
export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37
export X=~/rds/projects/covid/ace2/interval_genetic_data/interval_imputed_data/
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
export TMPDIR=${HPC_WORK}/work

function fastGWAsetup()
{
  cut -f1 ${analysis}/output/caprion${suffix}.pheno | sed '1d' > ${analysis}/output/caprion${suffix}.id
  paste -d' ' ${analysis}/output/caprion${suffix}.id ${analysis}/output/caprion${suffix}.id | \
  grep -v NA > ${analysis}/output/caprion${suffix}.id2
  seq 22 | \
  xargs -I {} echo ${analysis}/bgen/chr{}.bgen > ${analysis}/bgen/caprion.bgenlist
  echo ${analysis}/bgen/chr{1..22}.bgen | tr ' ' '\n' > ${analysis}/bgen/caprion${suffix}.bgenlist
  cat <(head -2 ${caprion}/data/merged_imputation.sample | awk '{if(NR==1) {print $0, "sex"} else {print $0, "D"}}') \
      <(grep -f ${analysis}/output/caprion${suffix}.id -w ${caprion}/data/merged_imputation.sample | \
        cut -d' ' -f1-2 | join - ${caprion}/data/merged_imputation.missing | \
        awk '{print $0, "NA"}') > ${analysis}/bgen/caprion.sample
  gcta-1.9 --bfile ${caprion}/data/merged_imputation --keep ${analysis}/output/caprion${suffix}.id2 \
           --make-grm --out ${analysis}/output/caprion${suffix} --threads 10
# a sparse GRM from SNP data
  gcta-1.9 --grm ${analysis}/output/caprion${suffix} --make-bK-sparse 0.05 --out ${analysis}/output/caprion-spgrm${suffix} --threads 10
# fastGWA mixed model
  cut -f1-38 --complement ${analysis}/output/caprion${suffix}.pheno | \
  head -1 | \
  tr '\t' '\n' > ${analysis}/output/caprion${suffix}.varlist
  if [ ! -d ~/Caprion/analysis/METAL${suffix}/qqmanhattanlz/slurm ]; then mkdir ~/Caprion/analysis/METAL${suffix}/qqmanhattanlz/slurm; fi
  if [ ! -d ~/Caprion/analysis/METAL${suffix}/pairs/slurm ]; then mkdir -p ~/Caprion/analysis/METAL${suffix}/pairs/slurm; fi
  if [ ! -d ~/Caprion/analysis/METAL${suffix}/miamiplot/slurm ]; then mkdir -p ~/Caprion/analysis/METAL${suffix}/miamiplot/slurm; fi
}
# slow and unnecessary
gcta-1.9 --mbgen ${analysis}/bgen/caprion.bgenlist --sample ${analysis}/bgen/caprion.sample \
         --make-grm --out ${analysis}/bgen/caprion --threads 10

# sbatch --export=ALL ${analysis}/5_pgwas.sb

function X()
{
  awk '!/NA/{print $1"_"$1}' ${analysis}/output/caprion${sufix}.id | \
  bcftools view -S - --force-samples ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -O z -o ${analysis}/output/INTERVAL-X.vcf.gz
  bcftools query -l ${analysis}/output/INTERVAL-X.vcf.gz | \
  tr '_' ' ' | \
  awk '{print $1}' | \
  bcftools reheader -s - ${analysis}/output/INTERVAL-X.vcf.gz -o ${analysis}/output/X.vcf.gz --threads 12
  bcftools index -tf ${analysis}/output/X.vcf.gz
  bcftools query -l ${analysis}/output/X.vcf.gz | awk '{print $1,$1}' > ${analysis}/output/chrX.idlist
  bcftools query -f "%ID\n" ${analysis}/output/X.vcf.gz > ${analysis}/output/chrX.snplist
  plink2 --vcf ${analysis}/output/X.vcf.gz --export bgen-1.2 bits=8 --double-id \
         --set-all-var-ids @:#_\$1_\$2 --new-id-max-allele-len 680 \
         --out ${analysis}/bgen/chrX
  bgenix -g ${analysis}/bgen/chrX.bgen -index -clobber
  plink2 --bgen ${analysis}/bgen/chrX.bgen ref-unknown --sample ${analysis}/bgen/chrX.sample \
         --freq \
         --out ${analysis}/bgen/chrX-freq
  sed '1d' ${analysis}/bgen/chrX-freq.afreq | \
  cut -f2 > ${analysis}/bgen/chrX.snplist
  cut -d' ' -f1 ${analysis}/output/caprion${suffix}-1.id | grep -f - ${analysis}/output/chrX.idlist \
     > ${analysis}/output/chrX${suffix}-1.id
  cut -d' ' -f1 ${analysis}/output/caprion${suffix}-2.id | grep -f - ${analysis}/output/chrX.idlist \
     > ${analysis}/output/chrX${suffix}-2.id
  cut -d' ' -f1 ${analysis}/output/caprion${suffix}-3.id | grep -f - ${analysis}/output/chrX.idlist \
     > ${analysis}/output/chrX${suffix}-3.id
  (
    parallel -C' ' 'cat ${analysis}/bgen/chr{}.snplist' ::: $(echo {1..22} X)
  )  | \
  uniq > ${analysis}/bgen/caprion.snplist
  awk 'NR>1 && $5>=0.01 && $5<=0.99 {print $2}' ${analysis}/bgen/chrX-freq.afreq > ${analysis}/bgen/chrX-0.01.snplist
  (
    parallel -C' ' 'cat ${analysis}/bgen/chr{}-0.01.snplist' ::: $(echo {1..22} X)
  ) | \
  uniq > ${analysis}/bgen/caprion-0.01.snplist
}

function lrlist()
{
  ls ${analysis}/work/*.fastGWA.gz | \
  xargs -l basename | \
  sed 's/.fastGWA.gz//g' | \
  grep -v chrX | \
  grep -v -f - ${analysis}/work/caprion${suffix}.varlist > ${analysis}/work/caprion${suffix}.lrlist
}

function collect()
{
  parallel -j10 --env caprion -C' ' '
    export caprion_protein={1}
    if [ -f ${analysis}/work/{2}-{1}.fastGWA.gz ] && [ -f ${analysis}/work/{2}-{1}-chrX.fastGWA.gz ]; then
       export out=${caprion_protein}.txt.bgz
       (
         echo -e "CHR\tSNP\tPOS\tEFF_ALLELE\tOTHER_ALLELE\tN\tEFF_ALLELE_FREQ\tBETA\tSE\tP"
         gunzip -c ${analysis}/work/spa/{2}-{1}.fastGWA.gz | \
         sed "1d"
         gunzip -c ${analysis}/work/spa/{2}-{1}-chrX.fastGWA.gz | \
         sed "1d"
       ) | \
       awk -vOFS="\t" "{print \$2,\$1,\$3,\$6,\$4,\$5,\$7,\$8,\$9,\$10}" | \
       bgzip -f > ${analysis}/work/${out}
       tabix -f -S1 -s2 -b3 -e3 ${analysis}/work/${out}
    fi
  ' ::: $(cat ${analysis}/work/caprion${suffix}.varlist)
  parallel -C' ' '
    export chr=chr{}
    awk "NR>1 && \$5>=0.001 && \$5<=0.999 {print \$2,\$3,\$4,\$5}" ${analysis}/bgen/${chr}-freq.afreq > ${analysis}/bgen/${chr}.freq
    awk "NR>1 && \$5>=0.01 && \$5<=0.99 {print \$2,\$3,\$4,\$5}" ${analysis}/bgen/${chr}-freq.afreq > ${analysis}/bgen/${chr}-0.01.freq
  ' ::: $(echo {1..22} X)
  cat ${analysis}/bgen/chr{1..22}.freq ${analysis}/bgen/chrX.freq > ${analysis}/bgen/caprion.freq
  cat ${analysis}/bgen/chr{1..22}-0.01.freq ${analysis}/bgen/chrX-0.01.freq > ${analysis}/bgen/caprion-0.01.freq
}

function tableMAF()
{
  seq 22 | \
  xargs -I {} cut -f15 ${ref}/impute_{}_interval.snpstats | \
  grep -v MAF | \
  sed '/^$/d' | \
  Rscript -e '
    maf <- scan("stdin");
    options(width=200)
    common_cutoffs <- cut(maf,c(0,0.01,0.05,0.5),include.lowest=TRUE);
    table(common_cutoffs);
    table(common_cutoffs)/length(maf);
    full_cutoffs <- cut(maf,c(0,0.0001,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5),include.lowest=TRUE);
    table(full_cutoffs);
    table(full_cutoffs)/length(maf)
  '
}

function EUR_1KGp3_chrX()
{
  export EUR=~/rds/public_databases/1000G
  export ldfile=${HPC_WORK}/locuszoom_1.4/data/1000G/genotypes/2014-10-14/EUR/chrX
  grep -w EUR ${EUR}/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 | \
  bcftools view -S - -o - -O z ${EUR}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz | \
  bcftools annotate --set-id 'chr%CHROM:%POS' -o ~/Caprion/analysis/work/chrX.vcf.gz -O z -
  plink --vcf ~/Caprion/analysis/work/chrX.vcf.gz --make-bed --out ${ldfile}
  sed -i 's/chrX:/chr23/' ${ldfile}.bim
}
# a bit clumsy with errors
# qctool -filetype vcf -g ~/Caprion/analysis/work/chrX.vcf.gz -ofiletype binary_ped -og ${ldfile}
# include.lowest=FALSE (default)
# Read 87696888 items
# common_cutoffs
#    (0,0.01] (0.01,0.05]  (0.05,0.5]
#    73182469     2990789     7038912
# common_cutoffs
#    (0,0.01] (0.01,0.05]  (0.05,0.5]
#  0.83449334  0.03410371  0.08026410
# full_cutoffs
#     (0,0.0001] (0.0001,0.001]   (0.001,0.01]    (0.01,0.05]     (0.05,0.1]      (0.1,0.2]      (0.2,0.3]      (0.3,0.4]      (0.4,0.5]
#       38883052       24922508        9376909        2990789        1414448        1825305        1393653        1236541        1168965
# full_cutoffs
#     (0,0.0001] (0.0001,0.001]   (0.001,0.01]    (0.01,0.05]     (0.05,0.1]      (0.1,0.2]      (0.2,0.3]      (0.3,0.4]      (0.4,0.5]
#     0.44338007     0.28418919     0.10692408     0.03410371     0.01612883     0.02081379     0.01589170     0.01410017     0.01332961
# include.lowest=TRUE
# Read 87696888 items
# common_cutoffs
#    [0,0.01] (0.01,0.05]  (0.05,0.5]
#    77667187     2990789     7038912
# common_cutoffs
#    [0,0.01] (0.01,0.05]  (0.05,0.5]
#  0.88563219  0.03410371  0.08026410
# full_cutoffs
#     [0,0.0001] (0.0001,0.001]   (0.001,0.01]    (0.01,0.05]     (0.05,0.1]      (0.1,0.2]      (0.2,0.3]      (0.3,0.4]      (0.4,0.5]
#       43367770       24922508        9376909        2990789        1414448        1825305        1393653        1236541        1168965
# full_cutoffs
#     [0,0.0001] (0.0001,0.001]   (0.001,0.01]    (0.01,0.05]     (0.05,0.1]      (0.1,0.2]      (0.2,0.3]      (0.3,0.4]      (0.4,0.5]
#     0.49451892     0.28418919     0.10692408     0.03410371     0.01612883     0.02081379     0.01589170     0.01410017     0.01332961

# by chromosome
# seq 22 | xargs -I {} bash -c "cut -f15 ${ref}/impute_{}_interval.snpstats > {}.MAF"

function uniprot_prot_list()
{
  Rscript -e '
    suppressMessages(library(dplyr))
    caprion <- select(pQTLtools::caprion,Protein,Accession,Gene) %>%
               mutate(Protein=gsub("_HUMAN","",Protein)) %>%
               rename(prot=Protein)
    write.table(caprion,file="~/Caprion/analysis/work/caprion.list",col.names=FALSE,row.names=FALSE,quote=FALSE)
 '
}

function caprion_ZYQ_classification()
{
  Rscript -e '
    library(dplyr)
    caprion_mc <- read.csv("~/Caprion/pilot/ZYQ_PC1_groups_20200703.csv") %>%
                  rename(PC1_group=pc1_group) %>%
                  mutate(PC1_group=as.factor(PC1_group))
    rownames(caprion_mc) <- with(caprion_mc,LIMS.ID)
    require(ggplot2)
    require(ggpubr)
    v <- ggplot(caprion_mc, aes(x=PC1_group, y=PC1, fill=PC1_group)) +
                geom_violin() +
                geom_boxplot(width=0.1) +
                xlab("PC1 group") +
                theme_minimal()
    means <- group_by(caprion_mc,PC1_group) %>%
             summarise(N=n(),Mean=signif(mean(PC1),3),SD=signif(sd(PC1),3))
    m <- ggtexttable(means, rows = NULL, theme = ttheme("mOrange"))
    p <- ggarrange(v,m,ncol=1,nrow=2,labels="2. ZYQ")
    ggsave("~/Caprion/analysis/work/ZYQ_grouping.png",width=10,height=12,units="in",device="png")
    suppressMessages(require(Biobase))
    suppressMessages(require(quantro))
    suppressMessages(require(sva))
    load("~/Caprion/pilot/ZYQ.rda")
  # Possible use
    g1 <- filter(caprion_mc,PC1_group=="group1")[["LIMS.ID"]]
    g2 <- filter(caprion_mc,PC1_group=="group2")[["LIMS.ID"]]
    g <- c(g1,g2)
    ZYQ_data1 <- protein_ZYQ[,g1]
    ZYQ_data2 <- protein_ZYQ[,g2]
    edata <- exprs(combine(ZYQ_data1,ZYQ_data2))
    png("~/Caprion/analysis/work/ZYQ_combat.png",res=300,width=20,height=20,units="in")
    combat_edata <- ComBat(dat=edata, batch=with(caprion_mc,PC1_group), par.prior=TRUE, ref.batch="group2", prior.plots=TRUE)
    dev.off()
    quantro_qtest  <- function(data,suffix,n.perm=10000,keep=TRUE)
    {
      qtest <- quantro(data,groupFactor=caprion_mc$PC1_group,B=n.perm)
      print(anova(qtest))
      print(quantroStat(qtest))
      qsp <- quantroStatPerm(qtest)
      print(quantroPvalPerm(qtest))
      if (keep) {
         png(paste0("~/Caprion/analysis/work/ZYQ_boxplot",suffix,".png"),res=300,width=50,height=12,units="in")
         matboxplot(data, groupFactor=with(caprion_mc,PC1_group), cex=0.4, pch=19, xaxt="n",
                    main="Boxplots of all proteins", xlab="Sample", ylab="Abundance level",
                    brewer.n=8, brewer.name="Dark2")
         legend('top', c("Group1", "Group2"), col = c(1, 2), lty = 1, lwd = 3)
         dev.off()
         save(qtest,file=paste0("~/Caprion/analysis/work/ZYQ_quantro",suffix,".rda"))
      }
    }
    quantro_qtest(edata,"")
    quantro_qtest(combat_edata,"_combat")
  '
}

function protein_mapping()
{
  cut -d, -f1-5 ~/Caprion/pre_qc_data/batch2/CAM1184-ZYQ/ZYQ_Comp_Neq1_Norm_Int_20200812.csv > ~/Caprion/pilot/data2/mapping_ZYQ.csv
  cut -f1-5 ~/Caprion/pre_qc_data/batch3/CAM1184-UDP/UDP_R1_Comp_Neq1_Raw_Int_Clean_20210412.txt > ~/Caprion/pilot/data3/mapping_UDP1.txt
  cut -f1-5 ~/Caprion/pre_qc_data/batch3/CAM1184-UDP/UDP_R2_Comp_Neq1_Raw_Int_Clean_20210412.txt > ~/Caprion/pilot/data3/mapping_UDP2.txt
  cut -f1-5 ~/Caprion/pre_qc_data/batch3/CAM1184-UDP/UDP_R3_Comp_Neq1_Raw_Int_Clean_20210412.txt > ~/Caprion/pilot/data3/mapping_UDP3.txt
  cut -d, -f1-5 ~/Caprion/pre_qc_data/batch4/UHZ_Comp_Raw_Int_20240118_v1.csv > ~/Caprion/pilot/data4/mapping_UHZ.csv
  Rscript -e '
    library(dplyr)
    library(VennDiagram)
    ZWK <- read.csv("~/Caprion/pre_qc_data/pilot/ZWK_protein_list.csv") %>% pull(Protein)
    ZYQ <- read.csv("~/Caprion/pre_qc_data/batch2/CAM1184-ZYQ/ZYQ_protein_list.csv") %>% pull(Protein)
    UDP <- read.csv("~/Caprion/pre_qc_data/batch3/CAM1184-UDP/UDP_protein_list.csv") %>% pull(Protein)
    load("~/Caprion/pilot/UHZ.rda")
    UHZ <- rownames(exprs(protein_UHZ))
    proteins <- list(ZWK,ZYQ,UDP,UHZ) %>%
                setNames(c("ZWK","ZYQ","UDP","UHZ"))
    venn.diagram(x=proteins,filename="~/Caprion/analysis/work/proteins.png",
                 disable.logging=TRUE, imagetype="png", output=TRUE,
                 height=12, width=12, units="cm", resolution=500,
                 fill=c("yellow","purple","green","blue"), otation.degree = 0)
    ZWK_ZYQ_UDP <- c(ZWK,ZYQ,UDP) %>% unique
    d1 <- setdiff(UHZ,ZWK_ZYQ_UDP) %>% sort %>% unique
    d2 <- setdiff(ZWK_ZYQ_UDP,UHZ) %>% sort %>% unique
    d1d2 <- cbind(d1,d2)
    d1d2[(length(d2)+1):length(d1),2] <- NA
    write.csv(d1d2,file="~/Caprion/analysis/work/diffs.csv")
)'
}
