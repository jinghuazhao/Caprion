#!/usr/bin/bash

function setup()
{
  for d in sentinels/slurm work qqmanhattanlz
  do
    if [ ! -d ${root}/${d} ]; then mkdir -p ${root}/${d}; fi
  done
}

function sb()
{
cat <<'EOL'> ${root}/${protein}-merge.sb
#!/usr/bin/bash

#SBATCH --job-name=_merge
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio

#SBATCH --export ALL
#SBATCH --array=1-_N_
#SBATCH --output=ROOT/sentinels/slurm/_merge_%A_%a.o
#SBATCH --error=ROOT/sentinels/slurm/_merge_%A_%a.e

export TMPDIR=${HPC_WORK}/work
export protein=PROTEIN
export isotope=$(head -1 ${root}/${protein}.pheno | awk -vn=${SLURM_ARRAY_TASK_ID} '{print $(n+2)}')

function pgz()
# 1. extract all significant SNPs
{
# zcat METAL/{}-1.tbl.gz | head -1
  zcat ${root}/METAL/${isotope}-1.tbl.gz | \
  awk 'NR>1 && $6>=0.01 && $12<=log(5e-8)/log(10)' | \
  sort -k1,1n -k2,2n | \
  gzip -f > ${root}/sentinels/${isotope}.p.gz
}

function _HLA()
# 2. handling HLA
{
  (
    zcat ${root}/METAL/${isotope}-1.tbl.gz | awk -vOFS="\t" 'NR==1{$1="Chrom";$2="Start" "\t" "End";print}'
    zcat ${root}/sentinels/${isotope}.p.gz | \
    awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
    awk '!($1 == "chr6" && $3 >= 25392021 && $3 < 33392022)'
    zcat ${root}/sentinels/${isotope}.p.gz | \
    awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
    awk '$1 == "chr6" && $3 >= 25392021 && $3 < 33392022' | \
      sort -k13,13g | \
      awk 'NR==1'
  ) > ${root}/sentinels/${isotope}_nold.p
  export lines=$(wc -l ${root}/sentinels/${isotope}_nold.p | cut -d' ' -f1)
  if [ $lines -eq 1 ]; then
     echo removing ${isotope}_nold with $lines lines
     rm ${root}/sentinels/${isotope}_nold.p
  fi
}

function sentinels()
{
  (
    mergeBed -i ${root}/sentinels/${isotope}_nold.p -d 1000000 -c 13 -o min | \
    awk -v OFS="\t" -v trait=${isotope} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "trait"
      print $0, trait
    }'
  ) > ${root}/sentinels/${isotope}.merged
  (
    cut -f1-4,13 ${root}/sentinels/${isotope}_nold.p| \
    bedtools intersect -a ${root}/sentinels/${isotope}.merged -b - -wa -wb | \
    awk '$4==$10' | \
    cut -f1-5,9,10 | \
    awk -v OFS="\t" -v trait=${isotope} '
    {
      if(NR==1) print "Chrom", "Start", "End", "P", "trait", "MarkerName", "CHR", "POS", "SNP", "P_check"
      chr=gsub(/chr/,"",$1)
      print $1,$2,$3,$4,trait,$1":"$3,chr,$3,$6,$7
    }'
  ) | uniq > ${root}/sentinels/${isotope}.sentinels

  Rscript -e '
    d <- file.path(Sys.getenv("root"),"sentinels")
    isotope <- Sys.getenv("isotope")
    f <- file.path(d,paste0(isotope,".sentinels"))
    m <- read.table(f,header=TRUE,as.is=TRUE)
    dim(m)
    head(m)
    suppressMessages(library(dplyr))
    t <- m %>% group_by(trait,Chrom,Start,End) %>% slice(which.min(P))
    t
    P <- with(m,P)
    p <- table(P)[table(P)>1]
    print(p)
    m <- subset(t,MarkerName!=".")
    cols <- c(1:5,9)
    write.table(m[,cols],file=file.path(d,paste0(isotope,".signals")),row.names=FALSE,quote=FALSE,sep="\t")
  '
}

function mean_by_genotype_gen_sample()
{
  read isotope chr bp pqtl < <(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]+1{gsub(/23/,"X",$2);print $1,$2,$3,$4}' ${root}/${protein}.merge)
  export isotope=${isotope}
  export chr=${chr}
  export bp=${bp}
  export pqtl=${pqtl}
  for batch in {1..3}
  do
    export batch=${batch}
    export out=${root}/means/${root}-${batch}-${pqtl}
    if [ ! -f ${out}.dat ]; then
       plink-2 --bgen ${pilot}/work/chr${chr}.bgen ref-unknown \
               --sample ${analysis}/work/caprion.sample \
               --chr ${chr} --from-bp ${bp} --to-bp ${bp} \
               --keep ${pilot}/caprion-${batch}.id \
               --pheno ${pilot}/work/caprion-${batch}.pheno --pheno-name ${prot} \
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
     protein <- Sys.getenv("protein")
     isotope <- Sys.getenv("isotope")
     pqtl <- Sys.getenv("pqtl")
     invisible(suppressMessages(sapply(c("dplyr","ggplot2","ggpubr"),require,character.only=TRUE)))
     process_batch <- function(batch, genotypes=c("100","010","001"))
     {
       datfile <- file.path(root,"means",paste(protein,batch,isotope,pqtl,sep="-"))
       dat <- read.table(paste0(datfile,".dat"),
                         colClasses=c("character","numeric","character","character","integer","character","character",rep("numeric",3)),
                         col.names=c("IID","Phenotype","chr","rsid","pos","A1","A2","g1","g2","g3")) %>%
              mutate(g=paste0(round(g1),round(g2),round(g3)),
                     Genotype=as.factor(case_when(g == genotypes[1] ~ paste0(A1,"/",A1),
                                                  g == genotypes[2] ~ paste0(A1,"/",A2),
                                                  g == genotypes[3] ~ paste0(A2,"/",A2),
                                                  TRUE ~ "---")))
       means <- group_by(dat,Genotype) %>%
                summarise(N=n(),Mean=mean(Phenotype))
       invisible(list(dat=dat,means=means))
     }
     v <- m <- list()
     for (batch in 1:3)
     {
         x <- process_batch(batch)
         v[[batch]] <- ggplot(with(x,dat), aes(x=Genotype, y=Phenotype, fill=Genotype)) +
                       geom_violin() +
                       geom_boxplot(width=0.1) +
                       theme_minimal()
         m[[batch]] <- ggtexttable(with(x,means), rows = NULL, theme = ttheme("mOrange"))
     }
     p <- ggarrange(v[[1]],v[[2]],v[[3]],m[[1]],m[[2]],m[[3]],ncol=3,nrow=2,labels=c("1. ZWK","2. ZYQ","3. UDP"))
     ggsave(file.path(root,"means",paste0(prot,"-",pqtl,".png")),device="png",width=16, height=10, units="in")
  '
}

function mean_by_genotype_dosage()
{
  read isotope chr bp pqtl < <(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]+1{gsub(/23/,"X",$2);print $1,$2,$3,$4}' ${analysis}/work/caprion.merge)
  export isotope=${isotope}
  export chr=${chr}
  export bp=${bp}
  export pqtl=${pqtl}
  for batch in {1..3}
  do
    export batch=${batch}
    export out=${root}/${protein}-${batch}-${isotope}-${pqtl}
    if [ ! -f ${out}.dat ]; then
       plink-2 --bgen ${pilot}/work/chr${chr}.bgen ref-unknown \
               --sample ${analysis}/work/caprion.sample \
               --chr ${chr} --from-bp ${bp} --to-bp ${bp} \
               --keep ${pilot}/work/caprion-${batch}.id \
               --pheno ${pilot}/work/caprion-${batch}.pheno --pheno-name ${prot} \
               --recode A include-alt \
               --out ${out}
       rm ${out}.log
    fi
  done
  Rscript -e '
     options(width=120)
     root <- Sys.getenv("root")
     protein <- Sys.getenv("protein")
     pqtl <- Sys.getenv("pqtl")
     invisible(suppressMessages(sapply(c("dplyr","ggplot2","ggpubr"),require,character.only=TRUE)))
     process_batch <- function(batch,digits=3)
     {
       datfile <- file.path(root,"means",paste("caprion",batch,prot,pqtl,sep="-"))
       dat <- read.delim(paste0(datfile,".raw"),check.names=FALSE,
                         colClasses=c("character","character","character","character","integer","numeric","numeric"))
       n7 <- names(dat)[7]
       names(dat)[6:7] <- c("Phenotype","Genotype")
       dat <- mutate(dat,Genotype=as.character(round(Genotype)))
       means <- group_by(dat,Genotype) %>%
                summarise(N=n(),Mean=signif(mean(Phenotype),digits))
       invisible(list(dat=dat,means=means,id=n7))
     }
     v <- m <- list()
     for (batch in 1:3)
     {
         x <- process_batch(batch)
         v[[batch]] <- ggplot(with(x,dat), aes(x=Genotype, y=Phenotype, fill=Genotype)) +
                       geom_violin() +
                       geom_boxplot(width=0.1) +
                       xlab(with(x,id)) +
                       theme_minimal()
         m[[batch]] <- ggtexttable(with(x,means), rows = NULL, theme = ttheme("mOrange"))
     }
     p <- ggarrange(v[[1]],v[[2]],v[[3]],m[[1]],m[[2]],m[[3]],ncol=3,nrow=2,labels=c("1. ZWK","2. ZYQ","3. UDP"))
     ggsave(file.path(root,"means",paste0(protein,"-",pqtl,".png")),device="png",width=16, height=10, units="in")
  '
}

# mean_by_genotype_dosage

for cmd in pgz _HLA sentinels; do $cmd; done
EOL

export N=$(head -1 ${root}/${protein}.pheno | awk '{print NF-2}')
sed -i "s|ROOT|${root}|;s|PROTEIN|${protein}|;s|_N_|${N}|" ${root}/${protein}-merge.sb

#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=cclake
}

function signals()
(
  cat ${root}/sentinels/*signals | \
  head -1 | \
  awk -v FS="\t" '{print "isotope",$0}'
  head -1 ${root}/${protein}.pheno | cut -d' ' -f1-2 --complement | tr ' ' '\n' | \
  parallel -C' ' '
    if [ -f ${root}/sentinels/{}.signals ]; then
       awk -v FS="\t" -v isotope={} "NR>1 {print isotope,\$0}" ${root}/sentinels/{}.signals
    fi
    if [ -f ${root}/sentinels/{}-chrX.signals ]; then
       awk -v FS="\t" -v isotope={} "NR>1 {print isotope,\$0}" ${root}/sentinels/{}-chrX.signals
    fi
  '
) > ${root}/${protein}.signals

function merge()
{
  cat <(gunzip -c ${root}/METAL/*-1.tbl.gz | head -1 | paste <(echo prot) -) \
      <(sed '1d;s/\t/ /g' ${root}/${protein}.signals | grep -v 'X:' | \
        parallel -C' ' -j20 'zgrep -w {7} ${root}/METAL/{1}-1.tbl.gz | paste <(echo {1}) -') \
      <(sed '1d;s/\t/ /g' ${root}/${protein}.signals | grep 'X:' | \
        parallel -C' ' -j20 'zgrep -w {7} ${root}/METAL/{1}-chrX-1.tbl.gz | paste <(echo {1}) -') \
      > ${root}/${protein}.merge
  cut -f11,14 ${root}/${protein}.merge | sed '1d' | awk -vOFS="\t" '{printf $2" "; if($1<0) print "-"; else print "+"}' | \
  sort -k1,1 -k2,2 | uniq -c | awk -vOFS="\t" '{print $1,$2,$3}'> ${root}/${protein}.dir
}

function cistrans()
{
  Rscript -e '
    options(width=120)
    suppressMessages(library(dplyr))
    suppressMessages(library(gap))
  # Directions
    root <- Sys.getenv("root")
    protein <- Sys.getenv("protein")
    caprion.dir <- within(read.table(paste0(root,"/",protein,".dir"),col.names=c("Count","Direction","Final")),{Direction=gsub(""," ",Direction)})
    knitr::kable(caprion.dir)
  # cis/trans classification
    signals <- read.table(paste0(root,"/",protein,".signals"),header=TRUE) %>% mutate(prot=protein)
    merged <- read.delim(paste0(root,"/",protein,".merge")) %>% mutate(prot=protein)
    names(merged)[1:4] <- c("prot","Chr","bp","SNP")
  # glist-hg19
    INF <- Sys.getenv("INF")
    glist_hg19 <- read.table(file.path(INF,"csd3","glist-hg19"),col.names=c("chr","start","end","gene")) %>%
                  filter(! chr %in% c("XY","Y"))
    X <- with(glist_hg19,chr=="X")
    glist_hg19[X,"chr"] <- "23"
    ucsc <- transmute(pQTLdata::hg19Tables,chr=gsub("chr","",X.chrom),start=chromStart,end=chromEnd,Gene=hgncSym) %>%
            select(Gene,chr,start,end)
    X <- with(ucsc,chr=="X")
    ucsc[X,"chr"] <- "23"
    caprion <- select(pQTLdata::caprion,Protein,Accession,Gene) %>%
               mutate(Protein=gsub("_HUMAN","",Protein)) %>%
               rename(prot=Protein)
    quadruple <- function(d,label) data.frame(Gene=label,chr=min(d$chr),start=min(d$start),end=max(d$end))
    caprion[c(55,237,433,435),]
    subset(glist_hg19,grepl("^AMY1",gene))
    subset(glist_hg19,grepl("^C4B",gene)&chr=="6")
    subset(glist_hg19,grepl("^HIST1H4|^HIST2H4[A-B]|HIST4H4",gene))
    subset(glist_hg19,grepl("^HBA",gene))
    AMY <- quadruple(subset(glist_hg19,grepl("^AMY1",gene)),label="AMY")
    C4B <- quadruple(subset(glist_hg19,grepl("^C4B",gene)&chr=="6"),label="C4B")
    HIST <- quadruple(subset(glist_hg19,grepl("^HIST1H4|^HIST2H4[A-B]|HIST4H4",gene)),label="HIST")
    HBA <- quadruple(subset(glist_hg19,grepl("^HBA",gene)),label="HBA")
    caprion_modified <- caprion
    caprion_modified[55,"Gene"] <- "AMY"
    caprion_modified[237,"Gene"] <- "C4B"
    caprion_modified[278,"Gene"] <- "C1orf123"
    caprion_modified[385,"Gene"] <- "FYB"
    caprion_modified[390,"Gene"] <- "FAM198A"
    caprion_modified[433,"Gene"] <- "HIST"
    caprion_modified[435,"Gene"] <- "HBA"
    caprion_modified[845,"Gene"] <- "SEPT2"
    APOC <- subset(glist_hg19,gene %in%c("APOC2","APOC4")) %>%
            rename(Gene=gene)
    ucsc_modified <- bind_rows(ucsc,APOC,AMY,C4B,HIST,HBA)
    pqtls <- select(merged,prot,SNP,log.P.) %>%
             mutate(log10p=-log.P.) %>%
             left_join(caprion_modified) %>%
             select(Gene,SNP,prot,log10p)
    posSNP <- select(merged,SNP,Chr,bp)
    cis.vs.trans <- qtlClassifier(pqtls,posSNP,ucsc_modified,1e6) %>%
                    mutate(geneChrom=as.integer(geneChrom),cis=if_else(Type=="cis",TRUE,FALSE))
    table(cis.vs.trans$Type)
    write.csv(cis.vs.trans,file=file.path(root,paste0(protein,".cis.vs.trans")),row.names=FALSE,quote=FALSE)
    png(file.path(root,"pqtl2d.png"),width=12,height=10,unit="in",res=300)
    r <- qtl2dplot(cis.vs.trans,chrlen=gap::hg19,snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                   gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                   TSS=TRUE,cis="cis",plot=TRUE,cex.labels=0.6,cex.points=0.6,
                   xlab="pQTL position",ylab="Gene position")
    dev.off()
    r <- qtl2dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                     snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                     gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                     TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                     xlab="pQTL position",ylab="Gene position")
    htmlwidgets::saveWidget(r,file=file.path(root,paste0(protein,"-pqtl2dplotly.html")))
    r <- qtl3dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                     snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                     gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                     TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                     xlab="pQTL position",ylab="Gene position")
    htmlwidgets::saveWidget(r,file=file.path(root,paste0(protein,"-pqtl3dplotly.html")))
  '
}

function fp()
{
  cp  ${root}/${protein}.merge ${root}/tbl.tsv
  cut -f2-4 ${root}/tbl.tsv | \
  awk 'NR>1' | \
  sort -k1,2n | \
  uniq | \
  awk -vOFS="\t" '{print $1":"$2,$3}' > ${root}/rsid.tsv
  (
    gunzip -c ${root}/${protein}-*fastGWA.gz | head -1
    awk 'NR>1' ${root}/tbl.tsv | \
    cut -f1,4,14 --output-delimiter=' ' | \
    parallel -j10 -C' ' '
      export direction=$(zgrep -w {2} ${root}/METAL/{1}-1.tbl.gz | cut -f13)
      let j=1
      for i in $(grep "Input File" ${root}/METAL/{1}-1.tbl.info | cut -d" " -f7)
      do
         export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
         if [ "$n" != "?" ]; then zgrep -H -w {2} $i; fi
         let j=$j+1
      done
  '
  ) | \
  sed 's/.gz//g' > ${root}/all.tsv
  Rscript -e '
    require(gap)
    require(dplyr)
    root <- Sys.getenv("root")
    protein <- Sys.getenv("protein")
    tbl <- read.delim(file.path(root,"tbl.tsv")) %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position)) %>%
           arrange(prot,SNP)
    all <- read.delim(file.path(root,"all.tsv")) %>%
           rename(EFFECT_ALLELE=A1,REFERENCE_ALLELE=A2) %>%
           mutate(CHR=gsub(paste0(root,"/",protein,"-|.fastGWA"),"",CHR)) %>%
           mutate(batch_prot_chr=strsplit(CHR,"-|:"),
                  batch=unlist(lapply(batch_prot_chr,"[",1)),
                  prot=unlist(lapply(batch_prot_chr,"[",2)),
                  CHR=unlist(lapply(batch_prot_chr,"[",3))) %>%
           mutate(MarkerName=paste0(CHR,":",POS),
                  study=case_when(batch == 1 ~ "1. ZWK",
                                  batch == 2 ~ "2. ZYQ",
                                  batch == 3 ~ "3. UDP",
                                  TRUE ~ "---")) %>%
           select(-batch_prot_chr)
    rsid <- read.table(file.path(root,"rsid.tsv"),col.names=c("MarkerName","rsid"))
    pdf(file.path(root,"fp.pdf"),width=10,height=8)
    METAL_forestplot(tbl,all,rsid)
    dev.off()
  '
}

function HetISq()
# Code extracted from caprion.Rmd
{
  Rscript -e '
    suppressMessages(require(dplyr))
    root <- Sys.getenv("root")
    protein <- Sys.getenv("protein")
    all <- read.delim(file.path(root,"all.tsv")) %>%
           mutate(CHR=gsub(paste0(root,"/",protein,"-|.fastGWA"),"",CHR)) %>%
           mutate(batch_pept_chr=strsplit(CHR,"-|:"),
                  batch=unlist(lapply(batch_pept_chr,"[",1)),
                  prot=unlist(lapply(batch_pept_chr,"[",2)),
                  CHR=unlist(lapply(batch_pept_chr,"[",3))) %>%
           mutate(MarkerName=paste0(CHR,":",POS),
                  Batch=case_when(batch == 1 ~ "1. ZWK",
                                  batch == 2 ~ "2. ZYQ",
                                  batch == 3 ~ "3. UDP",
                                  TRUE ~ "---"),
                  direction=case_when(sign(BETA) == -1 ~ "-", sign(BETA) == 1 ~ "+", sign(BETA) == 0 ~ "0", TRUE ~ "---")) %>%
           select(Batch,prot,-batch_pept_chr,MarkerName,SNP,A1,A2,N,AF1,BETA,SE,P,INFO,direction)
    b1 <- subset(all,Batch=="1. ZWK") %>% setNames(paste0(names(all),".ZWK")) %>% rename(prot=prot.ZWK, SNP=SNP.ZWK)
    b2 <- subset(all,Batch=="2. ZYQ") %>% setNames(paste0(names(all),".ZYQ"))
    b3 <- subset(all,Batch=="3. UDP") %>% setNames(paste0(names(all),".UDP"))
    b <- full_join(b1,b2,by=c('prot'='prot.ZYQ','SNP'='SNP.ZYQ')) %>%
         full_join(b3,by=c('prot'='prot.UDP','SNP'='SNP.UDP')) %>%
         mutate(directions=gsub("NA","?",paste0(direction.ZWK,direction.ZYQ,direction.UDP))) %>%
         select(-Batch.ZWK,-Batch.ZYQ,-Batch.UDP,direction.ZWK,direction.ZYQ,direction.UDP)
    tbl <- read.delim(file.path(root,"tbl.tsv")) %>%
           arrange(protein,MarkerName) %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position), index=1:n())
    Het <- filter(tbl,HetISq>=75) %>%
           mutate(prot=as.character(prot)) %>%
           select(prot,SNP,Direction,HetISq,index) %>%
           left_join(select(b,prot,SNP,P.ZWK,P.ZYQ,P.UDP,BETA.ZWK,BETA.ZYQ,BETA.UDP))
    write.csv(Het,file=file.path(root,"HetISq75.csv"),row.names=FALSE,quote=FALSE)
    write(Het[['index']],file=file.path(root,'HetISq75.index'),sep=",",ncolumns=nrow(Het))
  '
}

function qqmanhattan()
{
  module load python/3.7
  source ~/COVID-19/py37/bin/activate
  export dir=${root}/qqmanhattanlz
  head -1 ${root}/${protein}.pheno | cut -d' ' -f1-2 --complement | tr ' ' '\n' | \
  parallel -C' ' --env root '
  gunzip -c ${root}/METAL/{}-1.tbl.gz | \
  awk "{if (NR==1) print \"chromsome\",\"position\",\"log_pvalue\",\"beta\",\"se\"; else print \$1,\$2,-\$12,\$10,\$11}" > ${root}/work/{}.txt
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
       plot_title="{}" < ~/cambridge-ceu/turboman/turboman.r
  else
    cat <(echo chromosome position) \
        <(awk "NR>1{print \$1,\$2}" ${root}/sentinels/{}.signals) \
        > ${root}/work/{}.annotate
    R --slave --vanilla --args \
      input_data_path=${root}/work/{}.txt \
      output_data_rootname=${dir}/{}_manhattan \
      custom_peak_annotation_file_path=${root}/work/{}.annotate \
      reference_file_path=~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda \
      pvalue_sign=5e-8 \
      plot_title="{}" < ~/cambridge-ceu/turboman/turboman.r
    rm ${root}/work/{}.annotate
  fi
  rm ${root}/work/{}.txt
  convert +append ${dir}/{}_manhattan.png ${dir}/{}_qq.png -resize x500 -density 300 ${dir}/{}_qqmanhattan.png
  convert ${dir}/{}_qqmanhattan.png -quality 0 ${dir}/{}_qqmanhattan.jp2
  img2pdf -o ${dir}/{}_qqmanhattan.pdf ${dir}/{}_qqmanhattan.jp2
  rm ${dir}/{}_qqmanhattan.jp2
  '
  module load texlive
# pdfjam ${dir}/*_qqmanhattan.pdf --nup 1x1 --landscape --papersize '{7in,12in}' --outfile ${root}/qq+manhattan.pdf
  qpdf --empty --pages $(ls ${dir}/*_qqmanhattan.pdf) -- ${root}/qq+manhattan.pdf
  deactivate
}

function lz()
# for variant not in the reference panel, drop --refsnp, and then rename.
{
  module load python/2.7
  if [ -f ${root}/${protein}.signals ]; then
     (
       cat ${root}/${protein}.signals | \
       tr '\t' ' ' | \
       awk 'NR>1 && $2!="23" {print $6, $7}' | \
       parallel -j1 -C' ' --env root '
         zgrep -w {2} ${root}/METAL/{1}-1.tbl.gz | \
         awk -v isotope={1} -v rsid={2} "{print \$1,\$2-5e5,\$2+5e5,isotope,rsid}"
       '
     ) | \
     parallel -j1 -C ' ' --env root '
       (
         gunzip -c ${root}/METAL/{4}-1.tbl.gz | \
         awk -v OFS="\t" "NR==1 {\$12=\"log10P\";print}" | \
         cut -f1-3,12
         gunzip -c ${root}/METAL/{4}-1.tbl.gz | \
         awk -v chr={1} -v start={2} -v end={3} -v OFS="\t" "NR>1 && \$1==chr && \$2>=start && \$2<end {\$12=-\$12;print \$1,\$2,\$3,\$12}" | \
         sort -k1,1n -k2,2n
       ) > ${root}/work/{4}-{5}.lz
       locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${root}/work/{4}-{5}.lz \
                 --delim tab title="{4}-{5}" \
                 --markercol MarkerName --pvalcol log10P --no-transform --chr {1} --start {2} --end {3} --cache None \
                 --no-date --plotonly --prefix={4} --rundir ${root}/qqmanhattanlz --svg --refsnp {5}
       rm ${root}/work/{4}-{5}.lz
     '
  fi
  if [ -f ${root}/${protein}.signals ]; then
     (
       cat ${root}/${protein}.signals | \
       tr '\t' ' ' | \
       awk 'NR>1 && $2=="23" {print $6,$7}' ${root}/${protein}.signals | \
       parallel -j1 -C' ' --env root '
         zgrep -w {2} ${root}/METAL/{1}-chrX-1.tbl.gz | \
         awk -v isotope={1} -v rsid={2} "{print \$1,\$2-5e5,\$2+5e5,isotope,rsid}"
       '
     ) | \
     parallel -j1 -C ' ' --env root '
       (
         gunzip -c ${root}/METAL/{4}-chrX-1.tbl.gz | \
         awk -v OFS="\t" "NR==1 {\$12=\"log10P\";print}" | \
         cut -f1-3,12
         gunzip -c ${root}/METAL/{4}-chrX-1.tbl.gz | \
         awk -v chr={1} -v start={2} -v end={3} -v OFS="\t" "NR>1 && \$1==chr && \$2>=start && \$2<end {\$12=-\$12;print \$1,\$2,\$3,\$12}}" | \
         sort -k1,1n -k2,2n | \
         sed "s/_[A-Z]*_[A-Z]*//" | cut -f1-3,12 | sed "s/X/chr23/"
       ) > ${root}/work/{4}-chrX-{5}.lz
       locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${root}/work/{4}-chrX-{5}.lz \
                 --delim tab title="{4}-chrX-{5}" \
                 --markercol MarkerName --pvalcol log10P --no-transform --chr {1} --start {2} --end {3} --cache None \
                 --no-date --plotonly --prefix={4}-chrX --rundir ${root}/qqmanhattanlz \
                 --svg --refsnp $(echo chr{5} | sed "s/_[A-Z]*_[A-Z]*//")
       rm ${root}/work/{4}-chrX-{5}.lz
     '
  fi
}

function fplz()
{
  module load texlive
  export metal=${root}/METAL
# HSPB1_rs114800762 is missing as dug by the following code.
  join -a1 <(sed '1d' ${root}/${protein}.merge | awk '{print $1"_"$4}' | sort -k1,1 ) \
           <(ls ${root}/qqmanhattanlz/*pdf | xargs -l basename -s .pdf | awk '{print $1,NR}') | \
  awk 'NF<2' | \
  sed 's/_/ /' | \
  parallel -C' ' 'ls ${root}/qqmanhattanlz/{1}*pdf'
# forest/locuszoom left-right format
  ulimit -n
  ulimit -S -n 2048
  qpdf --empty --pages $(sed '1d' ${root}/${protein}.merge | sort -k1,1 -k4,4 | cut -f1,4 --output-delimiter=' ' | \
                         parallel -C' ' 'ls $(echo ${root}/qqmanhattanlz/{1}_{2}.pdf | sed "s/:/_/")') -- ${root}/lz2.pdf
  export npages=$(qpdf -show-npages ${root}/lz2.pdf)
  qpdf --pages . 1-$npages:odd -- ${root}/lz2.pdf ${root}/lz.pdf
# Split files, note the naming scheme
  pdfseparate ${root}/lz.pdf ${root}/work/temp-%04d-lz.pdf
  pdfseparate ${root}/fp.pdf ${root}/work/temp-%04d-fp.pdf
# left-right with very small file size
# Combine the final pdf
  pdfjam ${root}/work/temp-*-*.pdf --nup 2x1 --landscape --papersize '{7in,16in}' --outfile ${root}/fp+lz.pdf
  rm ${root}/work/temp*pdf
# qpdf fp+lz.pdf --pages . $(cat HetISq75.index) -- HetISq75.pdf
  qpdf ${root}/fp+lz.pdf --pages . $(sed '1d' ${root}/${protein}.merge | \
  sort -k1,1 -k4,4 | awk '$15>=75{printf " "NR}' | sed 's/ //;s/ /,/g') -- ${root}/HetISq75.pdf
}

export TMPDIR=${HPC_WORK}/work
export pilot=~/Caprion/pilot
export analysis=~/Caprion/analysis
for i in 754 # 319 # 490 # $(seq 987)
do
  export SLURM_ARRAY_TASK_ID=${i}
  export protein=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${pilot}/work/caprion.varlist)
  export root=${analysis}/peptide/${protein}
  setup
  sb
  sbatch --wait ${root}/${protein}-merge.sb
  signals
  merge
  cistrans
  fp
  HetISq
  qqmanhattan
  lz
  fplz
done