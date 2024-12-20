liftRegion <- function(x, flanking = 1e6)
{
  gr <- GenomicRanges::GRanges(seqnames = x$chr,
                               ranges = IRanges::IRanges(start = x$start - flanking, end = x$end + flanking))
  seqlevelsStyle(gr) <- "UCSC"
  f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
  chain <- rtracklayer::import.chain(f)
  gr38 <- rtracklayer::liftOver(gr, chain)
  lifted_chromosomes <- seqnames(gr38)
  lifted_start <- min(start(gr38))
  lifted_end <- max(end(gr38))
  chr <- gsub("chr","",colnames(table(lifted_chromosomes)))
  region <- paste0(chr, ":", lifted_start, "-", lifted_end)
  invisible(list(chr = chr, start = lifted_start, end = lifted_end, region = region))
}

sumstats <- function(prot,chr,region37,chain)
{
  cat("GWAS sumstats\n")
  tbl <- file.path(analysis,"METAL_dr",paste0(prot,"_dr-1.tbl.gz"))
  gwas_texts <- seqminer::tabix.read(tbl, tabixRange = region37)
  gwas_stats <- read.table(text = gwas_texts, sep = "\t", header = FALSE) %>%
                setNames(c("Chromosome","Position","ID","Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq",
                           "Effect","StdErr","logP","Direction","HetISq","HetChiSq","HetDf","logHetP","N"))
  gwas_granges <- with(gwas_stats,GRanges(seqnames = paste0("chr",dplyr::if_else(Chromosome==23,"X",paste(Chromosome))),
                       ranges = IRanges(start = Position, end = Position),
                       id = ID,REF=Allele2,ALT=Allele1,AF=Freq1,ES=Effect,SE=StdErr,LP=-logP,SS=N))
  f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
  chain <- rtracklayer::import.chain(f)
  gwas_stats_hg38 <- rtracklayer::liftOver(gwas_granges, chain) %>%
                     unlist() %>%
                     dplyr::as_tibble() %>%
                     dplyr::transmute(chromosome = seqnames,
                                      position = start, REF, ALT, AF, ES, SE, LP, SS) %>%
                     dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
                     dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
                     dplyr::group_by(id) %>%
                     dplyr::mutate(row_count = n()) %>%
                     dplyr::ungroup() %>%
                     dplyr::filter(row_count == 1) %>%
                     mutate(chromosome=gsub("chr","",chromosome))
}

microarray <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("a. eQTL datasets\n")
  microarray_df <- dplyr::filter(tabix_paths, quant_method == "microarray") %>%
                   dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(microarray_df$ftp_path), microarray_df$qtl_id[1])
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.CEDAR")
  column_names <- names(read.delim(hdr))
  summary_list <- purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region38,
                             selected_gene_id = ensGene, column_names))
  purrr::map_df(summary_list[lapply(summary_list,nrow)!=0],
                ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

rnaseq <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("b. Uniformly processed RNA-seq datasets\n")
  rnaseq_df <- dplyr::filter(tabix_paths, quant_method == "ge") %>%
               dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list, ~safe_import(., region38,
                             selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list[lapply(result_list,nrow)!=0], ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  purrr::map_df(result_list, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

gtex <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("c. GTEx_v8 imported eQTL datasets\n")
  fp <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_gtex.tsv")
  imported_tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>%
                          dplyr::mutate(ftp_path=file.path("~/rds/public_databases/GTEx/csv",basename(ftp_path)))
  gtex_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
             dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
  ftp_path_list <- setNames(as.list(gtex_df$ftp_path), gtex_df$qtl_id)
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list,
                             ~safe_import(., region38, selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list, ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0],
                                ~dplyr::filter(., !is.na(se)))
  prot <- sentinels[["prot"]][r]
  invisible(sapply(1:length(result_filtered), function(i) {
      if (!is.null(result_filtered[[i]])) {
          prot_name <- names(result_filtered)[i]
          if (!is.null(prot_name)) {
              f <- file.path(analysis, "coloc", "GTEx", "sumstats", paste0(prot, "-", prot_name, ".gz"))
              gz <- gzfile(f, "w")
              write.table(result_filtered[[i]], file = gz, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
              close(gz)
          } else {
              warning(paste("Missing name for index", i))
          }
      }
  }))
  purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
}

ge <- function(gwas_stats_hg38,ensGene,region38)
{
  cat("d. eQTL datasets\n")
  fp <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_ge.tsv")
  imported_tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>%
                          dplyr::mutate(ftp_path=file.path("~/rds/public_databases/eQTLCatalogue",basename(ftp_path)))
  ftp_path_list <- setNames(as.list(imported_tabix_paths$ftp_path), imported_tabix_paths$unique_id)
  hdr <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
  column_names <- names(read.delim(hdr))
  safe_import <- purrr::safely(import_eQTLCatalogue)
  summary_list <- purrr::map(ftp_path_list,
                             ~safe_import(., region38, selected_gene_id = ensGene, column_names))
  result_list <- purrr::map(summary_list, ~.$result)
  result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
  result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0],
                                ~dplyr::filter(., !is.na(se)))
  prot <- sentinels[["prot"]][r]
  invisible(sapply(1:length(result_filtered), function(i) {
      if (!is.null(result_filtered[[i]])) {
          prot_name <- names(result_filtered)[i]
          if (!is.null(prot_name)) {
              f <- file.path(analysis, "coloc", "eQTLCatalogue", "sumstats", paste0(prot, "-", prot_name, ".gz"))
              gz <- gzfile(f, "w")
              write.table(result_filtered[[i]], file = gz, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
              close(gz)
          } else {
              warning(paste("Missing name for index", i))
          }
      }
  }))
  purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "unique_id")
}

gtex_coloc <- function(prot,chr,ensGene,chain,region37,region38,out)
{
  gwas_stats_hg38 <- sumstats(prot,chr,region37,chain)
  df_gtex <- gtex(gwas_stats_hg38,ensGene,region38)
  if (!exists("df_gtex")) return
  saveRDS(df_gtex,file=paste0(out,".rds"))
  p <- ggplot(df_gtex, aes(x = PP.H4.abf)) + geom_histogram()
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggplot2::ggsave(plot = s, filename = paste0(out, ".assoc.pdf"), device = "pdf",
                  height = 15, width = 15, units = "cm", dpi = 300)
  ggplot2::ggsave(plot = p, filename = paste0(out, ".hist.pdf"), device = "pdf",
                  height = 15, width = 15, units = "cm", dpi = 300)
}

ge_coloc <- function(prot,chr,ensGene,chain,region37,region38,out)
{
  gwas_stats_hg38 <- sumstats(prot,chr,region37)
  df_ge <- ge(gwas_stats_hg38,ensGene,region38)
  if (!exists("df_ge")) return
  saveRDS(df_ge,file=paste0(out,".rds"))
  p <- ggplot(df_ge, aes(x = PP.H4.abf)) + geom_histogram()
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggsave(plot = s, filename = paste0(out, ".assoc.pdf"), device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
  ggsave(plot = p, filename = paste0(out, ".hist.pdf"), device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
}

all_coloc <- function(prot,chr,ensGene,chain,region37,region38,out)
{
  gwas_stats_hg38 <- sumstats(prot,chr,region37)
  df_microarray <- microarray(gwas_stats_hg38,ensGene,region38)
  df_rnaseq <- rnaseq(gwas_stats_hg38,ensGene,region38)
  df_gtex <- gtex(gwas_stats_hg38,ensGene,region38)
  df_ge <- ge(gwas_stats_hg38,ensGene,region38)
  if (exists("df_microarray") & exits("df_rnaseq") & exists("df_gtex") & exists("df_ge"))
  {
    coloc_df = dplyr::bind_rows(df_microarray, df_rnaseq, df_gtex, df_ge)
    saveRDS(coloc_df, file=paste0(out,"-all.rds"))
    p <- ggplot(coloc_df, aes(x = PP.H4.abf)) + geom_histogram()
  }
  s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
  ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
  ggsave(plot = p, filename = paste0(out, "-hist.pdf"), device = "pdf",
         height = 15, width = 15, units = "cm", dpi = 300)
}

single_run <- function(r, batch="GTEx")
{
  sentinel <- sentinels[r,]
  prot <- sentinel[["prot"]]
  chr <- sentinel[["geneChrom"]]
  ensRegion37 <- with(sentinel,
                      {
                        start <- geneStart-M
                        if (start<0) start <- 0
                        end <- geneEnd+M
                        paste0(chr,":",start,"-",end)
                      })
  ss <- subset(pQTLdata::caprion,Protein==paste0(sentinel[["prot"]],"_HUMAN"))
  ensGene <- ss[["ensGenes"]]
  x <- with(sentinel,list(chr=geneChrom,start=geneStart,end=geneEnd))
  lr <- liftRegion(x)
  ensRegion38 <- with(lr,paste0(chr,":",start-M,"-",end+M))
  cat(chr,ensGene,ensRegion37,ensRegion38,"\n")
  f <- file.path(analysis,"coloc",batch,with(sentinel,paste0(prot,"-",SNP)))
  if (batch=="GTEx")
  {
    gtex_coloc(prot,chr,ensGene,chain,ensRegion37,ensRegion38,f)
  } else {
    ge_coloc(prot,chr,ensGene,chain,ensRegion37,ensRegion38,f)
  }
}

collect <- function(batch="GTEx")
# to collect results when all single runs are done
{
  df_coloc <- data.frame()
  for(r in 1:nrow(sentinels))
  {
    prot <- sentinels[["prot"]][r]
    snpid <- sentinels[["SNP"]][r]
    rsid <- prot_rsid[["SNP"]][r]
    f <- file.path(analysis,"coloc",batch,paste0(prot,"-",snpid,".rds"))
    if (!file.exists(f)) next
    cat(prot,"-",rsid,"\n")
    rds <- readRDS(f)
    if (nrow(rds)==0) next
    df_coloc <- rbind(df_coloc,data.frame(prot=prot,rsid=rsid,snpid=snpid,rds))
  }
  f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
  chain <- rtracklayer::import.chain(f)
  gr <- with(arrange(pQTLdata::caprion,Protein), {
             GenomicRanges::GRanges(seqnames = chr,
                                    ranges = IRanges::IRanges(start = start - M, end = end + M),
                                    Protein=Protein,Gene=Gene)
             })
  seqlevelsStyle(gr) <- "UCSC"
  gr38 <- rtracklayer::liftOver(gr,chain)
  gr38_new <- lapply(gr38, function(gr) {
    valid_start <- start(gr)[!is.na(start(gr))]
    valid_end <- end(gr)[!is.na(end(gr))]
    if (length(valid_start) > 0 && length(valid_end) > 0) {
      min_start <- min(valid_start)
      min_end <- min(valid_end)
      new_gr <- GRanges(
        seqnames = seqnames(gr)[1],
        ranges = IRanges(start = min_start, end = min_end),
        strand = strand(gr)[1],
        Protein = ifelse(length(mcols(gr)$Protein) > 0, unique(mcols(gr)$Protein)[1], NA),
        Gene = ifelse(length(mcols(gr)$Gene) > 0, unique(mcols(gr)$Gene)[1], NA)
      )
      return(new_gr)
    } else {
      return(NULL)
    }
  })
  gr38_new <- gr38_new[!sapply(gr38_new, is.null)]
  df38 <- do.call(rbind, lapply(gr38_new, function(gr) {
             data.frame(
                seqnames = as.character(seqnames(gr)),
                start = start(gr),
                end = end(gr),
                strand = as.character(strand(gr)),
                Protein = mcols(gr)$Protein,
                Gene = mcols(gr)$Gene
             )
             })) %>%
          mutate(range38=paste0(gsub("chr","",seqnames),":",start,"-",end)) %>%
          select(Protein,Gene,range38)
  caprion_upd <- pQTLdata::caprion %>%
                 mutate(prot=gsub("_HUMAN","",Protein),gene=Gene) %>%
                 left_join(df38)
  df <- dplyr::rename(df_coloc,H0=PP.H0.abf,H1=PP.H1.abf,H2=PP.H2.abf,H3=PP.H3.abf,H4=PP.H4.abf) %>%
        dplyr::left_join(caprion_upd[c("prot","gene","range38")])
  if (batch=="GTEx") {
    df_coloc <- within(df,{qtl_id <- gsub("GTEx_V8_","",qtl_id)})
    write.table(subset(df,H4>=0.8),file=file.path(analysis,"coloc","GTEx.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    write.table(df,file=file.path(analysis,"coloc","GTEx-all.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    coloc <- merge(df_coloc,caprion_upd[c("prot","gene","range38")]) %>%
             mutate(prot,
                    H0=round(H0,2),
                    H1=round(H1,2),
                    H2=round(H2,2),
                    H3=round(H3,2),
                    H4=round(H4,2)) %>%
             setNames(c("Protein","Gene","RSid","SNPid","Tissue","nSNP","H0","H1","H2","H3","H4")) %>%
             select(Protein,Gene,RSid,Tissue,nSNP,H0,H1,H2,H3,H4)
    write.table(coloc,file=file.path(analysis,"coloc","GTEx-ST.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
  } else {
    write.table(subset(df,H4>=0.8),file=file.path(analysis,"coloc","eQTLCatalogue.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    write.table(df,file=file.path(analysis,"coloc","eQTLCatalogue-all.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
    eQTLCatalogue <- left_join(df,caprion_upd[c("prot","gene","range38")]) %>%
                     mutate(prot,
                            H0=round(H0,2),
                            H1=round(H1,2),
                            H2=round(H2,2),
                            H3=round(H3,2),
                            H4=round(H4,2)) %>%
                     setNames(c("Protein","Gene","RSid","SNPid","Study","nSNP","H0","H1","H2","H3","H4")) %>%
                     select(Protein,Gene,RSid,Study,nSNP,H0,H1,H2,H3,H4)
    write.table(eQTLCatalogue,file=file.path(analysis,"coloc","eQTLCatalogue-ST.tsv"),
                quote=FALSE,row.names=FALSE,sep="\t")
  }
}

loop_slowly <- function() for (r in 1:nrow(sentinels)) single_run(r)

# Environmental variables

pkgs <- c("dplyr", "gap", "ggplot2", "readr", "coloc", "GenomicRanges","pQTLtools","rtracklayer","seqminer")
invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))

options(width=200)
HOME <- Sys.getenv("HOME")
HPC_WORK <- Sys.getenv("HPC_WORK")
analysis <- Sys.getenv("analysis")
M <- 1e6
psum <- file.path(analysis,"coloc","sumstats")
if (!dir.exists(psum)) dir.create(psum)
gsum <- file.path(analysis,"coloc","GTEx","sumstats")
if (!dir.exists(gsum)) dir.create(gsum)
esum <- file.path(analysis,"coloc","eQTLCatalogue","sumstats")
if (!dir.exists(esum)) dir.create(esum)

sevens <- "
ENSG00000131142 - CCL25 19 8052318 8062660
ENSG00000125735 - TNFSF14 19 6661253 6670588
ENSG00000275302 - CCL4 17 36103827 36105621
ENSG00000274736 - CCL23 17 36013056 36017972
ENSG00000013725 - CD6 11 60971680 61020377
ENSG00000138675 - FGF5 4 80266639 80336680
ENSG00000277632 - CCL3 17 36088256 36090169
"
updates <- as.data.frame(scan(file=textConnection(sevens),what=list("","","",0,0,0))) %>%
           setNames(c("ensGenes","dash","gene","chromosome","start38","end38"))
caprion <- left_join(pQTLdata::caprion,updates)
sentinels <- subset(read.csv(file.path(analysis,"work","caprion_dr.cis.vs.trans")),cis)
f <- file.path(analysis,"work","snpid_dr.lst")
prot_rsid <- select(sentinels,prot,SNP) %>%
             dplyr::left_join(read.table(f,header=TRUE),by=c('SNP'='snpid')) %>%
             transmute(prot,SNP=dplyr::if_else(is.na(rsid)|rsid==".",SNP,rsid))
fp <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
tabix_paths <- read.delim(fp, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

r <- as.integer(Sys.getenv("r"))
single_run(r)
single_run(r,batch="eQTLCatalogue")

#collect()
#collect(batch="eQTLCatalogue")
