#!/usr/bin/bash

#SBATCH --job-name=_utils
#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --error=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/work/impute.e
#SBATCH --output=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/work/impute.o

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
export PERL5LIB=
module load ceuadmin/R
module load samtools/1.13/gcc/zwxn7ug3
module load perl/5.26.3_system/gcc-8.4.1-4cl2czq
module load libiconv/1.16/intel/64iicvbf
module load ceuadmin/ensembl-vep/111-icelake

export TMPDIR=${HPC_WORK}/work
export analysis=~/Caprion/analysis
export pre_qc_data=/rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics
export suffix=_dr

function impute()
{
   Rscript -e '
      suppressMessages(library(Biobase))
      suppressMessages(library(dplyr))
      suppressMessages(library(tibble))
      caprion <- "~/Caprion"
      load(file.path(caprion,"pilot","ZWK.rda"))
      load(file.path(caprion,"pilot","ZYQ.rda"))
      load(file.path(caprion,"pilot","UDP.rda"))
      load(file.path(caprion,"analysis","work","eSet.rda"))
      raw_ZYQ <- left_join(mapping_ZYQ,raw_ZYQ)
      raw_UDP <- left_join(mapping_UDP,raw_UDP)
      replace_below_threshold <- function(x, threshold = 50000) {
        x <- ifelse(x < threshold, NA, x)
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
      }
      normalize_minmax <- function(x) (x - min(x,na.rm=TRUE)) / (max(x,na.rm=TRUE) - min(x,na.rm=TRUE))
      normalize <- function(code)
      {
        dr <- Biobase::exprs(get(paste0("dr_",code))) %>% t()
        protein <- Biobase::exprs(get(paste0("protein_",code))) %>% t()
        proteins <- colnames(protein)
        peptide <- Biobase::exprs(get(paste0("peptide_",code))) %>% t()
        peptides <- colnames(peptide)
        raw <- get(paste0("raw_",code))
        samples <- grep(code,names(raw),value=TRUE)
        all_proteins <- unique(raw[["Protein"]])
        dropped_proteins <- grep("\\||-", all_proteins, value = TRUE)
      # 1,043 for ZYQ/UDP instead of 983/984
        isotopes <- subset(raw,Protein %in% setdiff(all_proteins,dropped_proteins)) %>%
                    dplyr::mutate(Isotope.Group.ID=as.character(Isotope.Group.ID))
        split_data <- split(isotopes,intersect(isotopes$Isotope.Group.ID,raw$Isotope.Group.ID))
        result_list <- sapply(isotopes[["Isotope.Group.ID"]], function(isotope) {
          pept_data <- split_data[[isotope]]
        # rownames(pept_data) <- isotope
          pept_data[samples] <- replace_below_threshold(unlist(pept_data[samples]))
          pept_data[paste0(samples, "_log2")] <- log2(unlist(pept_data[samples]+1))
          pept_data[paste0(samples, "_norm")] <- normalize_minmax(unlist(pept_data[samples]))
        # x <- pept_data[paste0(samples, "_norm")] |> unlist()
        # y <- pept_data[paste0(samples, "_log2")] |> unlist()
        # m <- lm(y~x)
        # pept_data[paste0(samples, "_log2s")] <- predict(m,newdata=data.frame(x),na.action=na.pass)
          return(pept_data)
        }, simplify = FALSE)
        result <- do.call(rbind, result_list)
        prot <- result[c("Protein",samples)] %>%
                dplyr::group_by(Protein) %>%
                dplyr::summarize(across(starts_with(code), sum, na.rm = TRUE), .groups = "drop") %>%
                tibble::column_to_rownames(var = "Protein") %>%
                t()
        z <-list(code=code,proteins=proteins,all_proteins=all_proteins,dropped_proteins=dropped_proteins,
                 peptide=peptide,peptides=peptides,dr=dr,protein=protein,proteins=proteins,
                 samples=samples,norm=result,prot=log2(prot+1))
        invisible(z)
      }
      impute_ZWK <- normalize("ZWK")
      save(impute_ZWK,file=file.path(caprion,"analysis","work","impute_ZWK.rda"))
      impute_ZYQ <- normalize("ZYQ")
      save(impute_ZYQ,file=file.path(caprion,"analysis","work","impute_ZYQ.rda"))
      impute_UDP <- normalize("UDP")
      save(impute_UDP,file=file.path(caprion,"analysis","work","impute_UDP.rda"))
   '
}

function mcpca()
{
   Rscript -e '
      suppressMessages(library(mclust))
      caprion <- "~/Caprion"
      for(batch in c("ZWK","ZYQ","UDP"))
      {
        load(file.path(caprion,"analysis","work",paste0("impute_",batch,".rda")))
        pdf(file.path(caprion,"analysis","work",paste0("impute_",batch,".pdf")))
        par(mfrow=c(2,2))
        attach(get(paste0("impute_",batch)))
        pca <- prcomp(protein,scale=TRUE)
        plot(pca$x[,1],pca$x[,2],xlab="PC1",ylab="PC2")
        pc1pc2 <- with(pca,x)[,1:2]
        mc <- Mclust(pc1pc2,G=4)
        scatterplot3d::scatterplot3d(with(pca,x[,c(1,2,3)]),color=c("blue","red")[mc$classification],
                                     main="Plot of the PC1, PC2 and PC3", pch=16)
        legend("right", legend=levels(as.factor(mc$classification)), col=c("blue", "red", "black","yellow"), pch=16)
        pca2 <- prcomp(prot[,proteins],scale=FALSE)
        plot(pca2$x[,1],pca2$x[,2],xlab="PC1",ylab="PC2")
        pc1pc2 <- with(pca2,x)[,1:2]
        mc <- Mclust(pc1pc2,G=2)
        scatterplot3d::scatterplot3d(with(pca2,x[,c(1,2,3)]),color=c("blue","red")[mc$classification],
                                     main="Plot of the PC1, PC2 and PC3", pch=16)
        legend("right", legend=levels(as.factor(mc$classification)), col=c("blue", "red"), pch=16)
        detach(get(paste0("impute_",batch)))
        dev.off()
      }
   '
}

