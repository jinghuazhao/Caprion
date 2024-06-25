#!/usr/bin/bash

function step2_pqtl_collect()
{
cat <<'EOL'> ${root}/${protein}-step2.sb
#!/usr/bin/bash

#SBATCH --job-name=_2-LABEL
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account PETERS-SL3-CPU
#SBATCH --partition icelake-himem

#SBATCH --export ALL
#SBATCH --output=ROOT/sentinels/slurm/_step2_LABEL.o
#SBATCH --error=ROOT/sentinels/slurm/_step2_LABEL.e

export TMPDIR=${HPC_WORK}/work
export protein=PROTEIN

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load ceuadmin/R
module load samtools/1.13/gcc/zwxn7ug3
module load perl/5.26.3_system/gcc-8.4.1-4cl2czq
module load libiconv/1.16/intel/64iicvbf
module load ceuadmin/ensembl-vep/111-icelake

function pqtls()
(
  cat ${root}/sentinels/*signals | \
  head -1 | \
  awk -v FS="\t" '{print "isotope",$0}'
  head -1 ${root}/${protein}.pheno | cut -d' ' -f1-2 --complement | tr ' ' '\n' | \
  parallel -C' ' '
    if [ -f ${root}/sentinels/{}.signals ]; then
       awk -v FS="\t" -v OFS="\t" -v isotope={} "NR>1 {print isotope,\$0}" ${root}/sentinels/{}.signals
    fi
  '
) > ${root}/${protein}.signals

function merge()
{
  cat <(gunzip -c ${root}/METAL/*-1.tbl.gz | head -1 | paste <(echo prot) -) \
      <(sed '1d;s/\t/ /g' ${root}/${protein}.signals | \
        parallel -C' ' -j10 'zgrep -w {7} ${root}/METAL/{1}-1.tbl.gz | paste <(echo {1}) -') \
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
    caprion.dir <- within(read.table(paste0(root,"/",protein,".dir"),col.names=c("Count","Direction","Final"),fill=TRUE),
                          {Direction=gsub(""," ",Direction)})
    knitr::kable(caprion.dir)
  # cis/trans classification
    signals <- read.table(paste0(root,"/",protein,".signals"),header=TRUE) %>% mutate(prot=protein)
    merged <- read.delim(paste0(root,"/",protein,".merge")) %>% mutate(isotope=prot,prot=protein)
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
    pqtls <- select(merged,prot,SNP,log.P.,isotope) %>%
             mutate(log10p=-log.P.) %>%
             left_join(caprion_modified) %>%
             filter(complete.cases(.)) %>%
             select(Gene,SNP,prot,log10p,isotope)
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

function vep_annotate()
{
  export cvt=${root}/${protein}.cis.vs.trans
  awk -v FS="," 'NR>1{print $5}' ${cvt} | \
  sort -k1,1 | \
  uniq | \
  parallel -C' ' '
    export isotope={}
    (
      echo "##fileformat=VCFv4.0"
      echo "#CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO"
      awk -vFS="," "\$5==ENVIRON[\"isotope\"] {print \$2}" ${cvt} | \
      sort -k1,1 | \
      zgrep -f - -w ${root}/METAL/{}-1.tbl.gz | \
      cut -f1-5 | \
      awk "{gsub(/23/,\"X\",\$1);print \$1,\$2,\$3,toupper(\$4),toupper(\$5),\".\",\".\",\".\"}"
    ) | \
    tr " " "\t" > ${root}/METAL/vep/{}.vcf
  # VEP annotation
    cd ${HPC_WORK}/loftee
    vep --input_file ${root}/METAL/vep/{}.vcf \
        --output_file ${root}/METAL/vep/{}.tab --force_overwrite \
        --cache --dir_cache /usr/local/Cluster-Apps/ceuadmin/ensembl-vep/111-icelake/.vep \
        --offline \
        --species homo_sapiens --assembly GRCh37 --pick --nearest symbol --symbol \
        --tab
    cd -
    (
      echo chromosome position nearest_gene_name cistrans
      awk -vFS="," "\$5==ENVIRON[\"isotope\"] {print \$2,\$9,\$10,\$11}" ${cvt} | \
      sort -k1,1 | \
      join - <(awk "!/#/{print \$1,\$21}" ${root}/METAL/vep/{}.tab | sort -k1,1) | \
      awk "{print \$2,\$3,\$5,\$4}" | \
      sort -k1,1n -k2,2n | \
      uniq
    ) > ${root}/METAL/vep/{}.txt
  '
}

for cmd in pqtls merge cistrans vep_annotate; do $cmd; done
EOL

sed -i "s|ROOT|${root}|;s|LABEL|${protein}|;s|PROTEIN|${protein}|" ${root}/${protein}-step2.sb
}

source 0_setup.sh

# all proteins:
while IFS=":" read -r protein_index protein; do
    export protein_index
    export protein
    echo ${protein_index} ${protein}
    export root=~/Caprion/analysis/peptide/${protein}
    export pheno=${analysis}/peptide/${protein}/${protein}.pheno
    export N=$(awk 'NR==1{print NF-2}' ${pheno})
  # 2. Collection of signals
    echo Step 2:
    step2_pqtl_collect
  # fplz should be here
    sbatch ${root}/${protein}-step2.sb
done < <(xargs -n 2 < ${analysis}/peptide_progs/benchmark2.names | \
         grep -n -f ${analysis}/peptide_progs/benchmark2.names -v -w ${varlist})

