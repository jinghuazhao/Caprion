#!/usr/bin/bash

function setup()
{
  if [ ! -d ${root}/sentinels/slurm ]; then mkdir -p ${root}/sentinels/slurm; fi
}

function merge_sb()
{
cat <<'EOL'> ${root}/merge.sb
#!/usr/bin/bash

#SBATCH --job-name=_merge
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio

#SBATCH --export ALL
#SBATCH --output=ROOT/sentinels/slurm/merge.o
#SBATCH --error=ROOT/sentinels/slurm/merge.e

export TMPDIR=${HPC_WORK}/work

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
    awk -v OFS="\t" -v trait=${p} '
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
  export caprion=~/Caprion
  read prot chr bp pqtl < <(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]+1{gsub(/23/,"X",$2);print $1,$2,$3,$4}' ${analysis}/peptide/caprion.merge)
  export prot=${prot}
  export chr=${chr}
  export bp=${bp}
  export pqtl=${pqtl}
  for batch in {1..3}
  do
    export batch=${batch}
    export out=${root}/means/${root}-${batch}-${pqtl}
    if [ ! -f ${out}.dat ]; then
       plink-2 --bgen ${caprion}/pilot/work/chr${chr}.bgen ref-unknown \
               --sample ${caprion}/analysis/work/caprion.sample \
               --chr ${chr} --from-bp ${bp} --to-bp ${bp} \
               --keep ${caprion}/pilot/work/caprion-${batch}.id \
               --pheno ${caprion}/pilot/work/caprion-${batch}.pheno --pheno-name ${prot} \
               --recode oxford \
               --out ${out}
       paste <(awk 'NR>2{print $1,$5}' ${out}.sample) \
             <(awk '{for(i=0;i<(NF-5)/3;i++) print $1,$2,$3,$4,$5, $(6+i),$(7+i),$(8+i)}' ${out}.gen) > ${out}.dat
       rm ${out}.gen ${out}.sample ${out}.log
    fi
  done
  Rscript -e '
     options(width=120)
     caprion <- Sys.getenv("caprion")
     prot <- Sys.getenv("prot")
     pqtl <- Sys.getenv("pqtl")
     invisible(suppressMessages(sapply(c("dplyr","ggplot2","ggpubr"),require,character.only=TRUE)))
     process_batch <- function(batch, genotypes=c("100","010","001"))
     {
       datfile <- file.path(caprion,"analysis","pgwas","means",paste("caprion",batch,prot,pqtl,sep="-"))
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
     ggsave(file.path(caprion,"analysis","pgwas","means",paste0(prot,"-",pqtl,".png")),device="png",width=16, height=10, units="in")
  '
}

function mean_by_genotype_dosage()
{
  export caprion=~/Caprion
  read prot chr bp pqtl < <(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]+1{gsub(/23/,"X",$2);print $1,$2,$3,$4}' ${caprion}/analysis/work/caprion.merge)
  export prot=${prot}
  export chr=${chr}
  export bp=${bp}
  export pqtl=${pqtl}
  for batch in {1..3}
  do
    export batch=${batch}
    export out=${caprion}/analysis/pgwas/means/caprion-${batch}-${prot}-${pqtl}
    if [ ! -f ${out}.dat ]; then
       plink-2 --bgen ${caprion}/pilot/work/chr${chr}.bgen ref-unknown \
               --sample ${caprion}/analysis/work/caprion.sample \
               --chr ${chr} --from-bp ${bp} --to-bp ${bp} \
               --keep ${caprion}/pilot/work/caprion-${batch}.id \
               --pheno ${caprion}/pilot/work/caprion-${batch}.pheno --pheno-name ${prot} \
               --recode A include-alt \
               --out ${out}
       rm ${out}.log
    fi
  done
  Rscript -e '
     options(width=120)
     caprion <- Sys.getenv("caprion")
     prot <- Sys.getenv("prot")
     pqtl <- Sys.getenv("pqtl")
     invisible(suppressMessages(sapply(c("dplyr","ggplot2","ggpubr"),require,character.only=TRUE)))
     process_batch <- function(batch,digits=3)
     {
       datfile <- file.path(caprion,"analysis","pgwas","means",paste("caprion",batch,prot,pqtl,sep="-"))
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
     ggsave(file.path(caprion,"analysis","pgwas","means",paste0(prot,"-",pqtl,".png")),device="png",width=16, height=10, units="in")
  '
}

# mean_by_genotype_dosage

for isotope in $(head -1 ${root}/${protein}.pheno | cut -d' ' -f1-2 --complement)
do
    for cmd in pgz _HLA sentinels; do $cmd; done
done
EOL

sed -i "s|ROOT|${root}|" ${root}/merge.sb

#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=cclake
}

function signals()
(
  cat ${root}/sentinels/*signals | \
  head -1 | \
  awk -v FS="\t" '{print "prot",$0}'
  head -1 ${root}/${protein}.pheno | cut -d' ' -f1-2 --complement | tr ' ' '\n' | \
  parallel -C' ' '
    if [ -f ${root}/sentinels/{}.signals ]; then
       awk -v FS="\t" -v prot={} "NR>1 {print prot,\$0}" ${root}/sentinels/{}.signals
    fi
    if [ -f ${root}/sentinels/{}-chrX.signals ]; then
       awk -v FS="\t" -v prot={} "NR>1 {print prot,\$0}" ${root}/sentinels/{}-chrX.signals
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
    signals <- read.table(paste0(root,"/",protein,".signals"),header=TRUE)
    merged <- read.delim(paste0(root,"/",protein,".merge"))
    names(merged)[1:4] <- c("prot","Chr","bp","SNP")
  # glist-hg19
    INF <- Sys.getenv("INF")
    glist_hg19 <- read.table(file.path(INF,"csd3","glist-hg19"),col.names=c("chr","start","end","gene")) %>%
                  filter(! chr %in% c("XY","Y"))
    X <- with(glist_hg19,chr=="X")
    glist_hg19[X,"chr"] <- "23"
    ucsc <- transmute(pQTLtools::hg19Tables,chr=gsub("chr","",X.chrom),start=chromStart,end=chromEnd,Gene=hgncSym) %>%
            select(Gene,chr,start,end)
    X <- with(ucsc,chr=="X")
    ucsc[X,"chr"] <- "23"
    caprion <- select(pQTLtools::caprion,Protein,Accession,Gene) %>%
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
    write.csv(cis.vs.trans,file=file.path(work,"caprion.cis.vs.trans"),row.names=FALSE,quote=FALSE)
    png(file.path(work,"caprion.pqtl2d.png"),width=12,height=10,unit="in",res=300)
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
    htmlwidgets::saveWidget(r,file=file.path(work,"caprion.pqtl2dplotly.html"))
    r <- qtl3dplotly(cis.vs.trans,chrlen=gap::hg19,qtl.id="SNP",qtl.prefix="pQTL:",target.type="Protein",
                     snp_name="SNP",snp_chr="SNPChrom",snp_pos="SNPPos",
                     gene_chr="geneChrom",gene_start="geneStart",gene_end="geneEnd",trait="prot",gene="Gene",
                     TSS=FALSE,cis="cis",cex.labels=0.6,cex.points=0.6,
                     xlab="pQTL position",ylab="Gene position")
    htmlwidgets::saveWidget(r,file=file.path(work,"caprion.pqtl3dplotly.html"))
  '
}

function mean()
{
  export caprion=~/Caprion
  awk '{gsub(/NA/,"0",$NF);print}' ${root}/caprion.sample > ${root}/caprion.sample
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
    gunzip -c ${root}/${root}-*fastGWA.gz | head -1
    awk 'NR>1' ${root}/tbl.tsv | \
    cut -f1,4,14 --output-delimiter=' ' | \
    parallel -j10 -C' ' '
      export direction=$(zgrep -w {2} ${analysis}/METAL/{1}-1.tbl.gz | cut -f13)
      let j=1
      for i in $(grep "Input File" ${analysis}/METAL/{1}-1.tbl.info | cut -d" " -f7)
      do
         export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
         if [ "$n" != "?" ]; then zgrep -H -w {2} $i; fi
         let j=$j+1
      done
  '
  ) | \
  sed 's|/data/jinhua/INF/sumstats||g;s/.gz//g' > ${root}/all.tsv
  Rscript -e '
    require(gap)
    require(dplyr)
    root <- Sys.getenv("root")
    tbl <- read.delim(file.path(root,"tbl.tsv")) %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position)) %>%
           arrange(prot,SNP)
    all <- read.delim(file.path(root,"all.tsv")) %>%
           rename(EFFECT_ALLELE=A1,REFERENCE_ALLELE=A2) %>%
           mutate(CHR=gsub("/home/jhz22/Caprion/analysis/work/pgwas/caprion-|.fastGWA","",CHR)) %>%
           mutate(batch_prot_chr=strsplit(CHR,"-|:"),
                  batch=unlist(lapply(batch_prot_chr,"[",1)),
                  prot=unlist(lapply(batch_prot_chr,"[",2)),
                  CHR=unlist(lapply(batch_prot_chr,"[",3))) %>%
           mutate(MarkerName=paste0(CHR,":",POS),
                  study=case_when(batch == batch[1] ~ "1. ZWK",
                                  batch == batch[2] ~ "2. ZYQ",
                                  batch == batch[3] ~ "3. UDP",
                                  TRUE ~ "---")) %>%
           select(-batch_prot_chr)
    rsid <- read.table("~/Caprion/analysis/work/rsid.tsv",col.names=c("MarkerName","rsid"))
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
    all <- read.delim(file.path(root,"all.tsv")) %>%
           mutate(CHR=gsub("/home/jhz22/Caprion/analysis/work/pgwas/caprion-|.fastGWA","",CHR)) %>%
           mutate(batch_prot_chr=strsplit(CHR,"-|:"),
                  batch=unlist(lapply(batch_prot_chr,"[",1)),
                  prot=unlist(lapply(batch_prot_chr,"[",2)),
                  CHR=unlist(lapply(batch_prot_chr,"[",3))) %>%
           mutate(MarkerName=paste0(CHR,":",POS),
                  Batch=case_when(batch == batch[1] ~ "1. ZWK",
                                  batch == batch[2] ~ "2. ZYQ",
                                  batch == batch[3] ~ "3. UDP",
                                  TRUE ~ "---"),
                  direction=case_when(sign(BETA) == -1 ~ "-", sign(BETA) == 1 ~ "+", sign(BETA) == 0 ~ "0", TRUE ~ "---")) %>%
           select(Batch,prot,-batch_prot_chr,MarkerName,SNP,A1,A2,N,AF1,BETA,SE,P,INFO,direction)
    b1 <- subset(all,Batch=="1. ZWK")
    names(b1) <- paste0(names(all),".ZWK")
    b1 <- rename(b1, prot=prot.ZWK, SNP=SNP.ZWK)
    b2 <- subset(all,Batch=="2. ZYQ")
    names(b2) <- paste0(names(all),".ZYQ")
    b3 <- subset(all,Batch=="3. UDP")
    names(b3) <- paste0(names(all),".UDP")
    b <- full_join(b1,b2,by=c('prot'='prot.ZYQ','SNP'='SNP.ZYQ')) %>% full_join(b3,by=c('prot'='prot.UDP','SNP'='SNP.UDP')) %>%
         mutate(directions=gsub("NA","?",paste0(direction.ZWK,direction.ZYQ,direction.UDP))) %>%
         select(-Batch.ZWK,-Batch.ZYQ,-Batch.UDP,direction.ZWK,direction.ZYQ,direction.UDP)
    tbl <- read.delim(file.path(root,"tbl.tsv")) %>%
           arrange(prot,MarkerName) %>%
           mutate(SNP=MarkerName,MarkerName=paste0(Chromosome,":",Position), index=1:n())
    Het <- filter(tbl,HetISq>=75) %>%
           select(prot,SNP,Direction,HetISq,index) %>%
           left_join(select(b,prot,SNP,P.ZWK,P.ZYQ,P.UDP,BETA.ZWK,BETA.ZYQ,BETA.UDP))
    write.csv(Het,file=file.path(root,"HetISq75.csv"),row.names=FALSE,quote=FALSE)
    write(Het[['index']],file=file.path(root,'HetISq75.index'),sep=",",ncolumns=nrow(Het))
  '
}

function fplz()
{
  export metal=~/Caprion/analysis/METAL
# HSPB1_rs114800762 is missing as dug by the following code.
  join -a1 <(sed '1d' work/caprion.merge | awk '{print $1"_"$4}' | sort -k1,1 ) \
           <(ls METAL/qqmanhattanlz/lz/*pdf | xargs -l basename -s .pdf | awk '{print $1,NR}') | \
  awk 'NF<2' | \
  sed 's/_/ /' | \
  parallel -C' ' 'ls METAL/qqmanhattanlz/lz/{1}*pdf'
# forest/locuszoom left-right format
  ulimit -n
  ulimit -S -n 2048
  qpdf --empty --pages $(sed '1d' ~/Caprion/analysis/work/caprion.merge | sort -k1,1 -k4,4 | cut -f1,4 --output-delimiter=' ' | \
                         parallel -C' ' 'ls $(echo METAL/qqmanhattanlz/lz/{1}_{2}.pdf | sed "s/:/_/")') -- lz2.pdf
  export npages=$(qpdf -show-npages lz2.pdf)
  qpdf --pages . 1-$npages:odd -- lz2.pdf lz.pdf
# Split files, note the naming scheme
  pdfseparate lz.pdf temp-%04d-lz.pdf
  pdfseparate ${metal}/fp/fp.pdf temp-%04d-fp.pdf
# left-right with very small file size
# Combine the final pdf
  pdfjam temp-*-*.pdf --nup 2x1 --landscape --papersize '{7in,16in}' --outfile fp+lz.pdf
  rm temp*pdf
# qpdf fp+lz.pdf --pages . $(cat HetISq75.index) -- HetISq75.pdf
  qpdf fp+lz.pdf --pages . $(sed '1d' caprion.merge | sort -k1,1 -k4,4 | awk '$15>=75{printf " "NR}' | sed 's/ //;s/ /,/g') -- HetISq75.pdf
}

for i in 14 # $(seq 987)
do
  export SLURM_ARRAY_TASK_ID=${i}
  export TMPDIR=${HPC_WORK}/work
  export pilot=~/Caprion/pilot
  export analysis=~/Caprion/analysis
  export protein=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${pilot}/work/caprion.varlist)
  export root=${analysis}/peptide/${protein}

  setup
  merge_sb
  sbatch --wait ${root}/merge.sb
  signals
  merge
  cistrans
  fp
done
