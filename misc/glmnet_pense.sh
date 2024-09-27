#!/usr/bin/bash

#SBATCH --job-name=_glmnet
#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --error=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/work/glmnet.e
#SBATCH --output=/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/work/glmnet.o

function glmnet_pense()
{
  Rscript -e '
    options(width=200)
    suppressMessages(library(Biobase))
    suppressMessages(library(glmnet))
    suppressMessages(library(pense))
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyr))
    suppressMessages(library(parallel))
    cluster <- makeCluster(5)
    analysis <- "~/Caprion/analysis"
    for (p in c("A1BG","APOB","EPCR","ERAP2","PROC"))
    {
       isotope <- read.table(file.path(analysis,"peptide",p,paste0(p,".pheno")),
                             check.names=FALSE,header=TRUE, nrows=1) %>%
                  select(-FID,-IID) %>%
                  names()
       pdf(file.path(analysis,"work",paste0(p,"-enet.pdf")),height=10,width=12)
       opar <- par()
       par(mfrow=c(2,1))
       isotope_ZWK <- isotope_ZYQ <- isotope_UDP <- isotope
       batches <- c("ZWK","ZYQ","UDP")
       sink(file.path(analysis,"work",paste0(p,"-enet.md")))
       for (batch in batches) {
          tryCatch({
            load(paste0("~/Caprion/pilot/",batch,".rda"))
            uniprot <- paste0(p,"_HUMAN")
            protein_batch <- get(paste0("protein","_",batch))
            peptide_batch <- get(paste0("peptide","_",batch))
            dr_batch <- get(paste0("dr","_",batch))
            mapping_batch <- get(paste0("mapping","_",batch))
            protein <- exprs(protein_batch)[rownames(exprs(protein_batch))==uniprot,]
            dr <- exprs(dr_batch)[rownames(exprs(dr_batch))==uniprot,]
            # isotope <- mapping_batch[mapping_batch$Protein == uniprot, "Isotope.Group.ID"]
            isotope <- get(paste0("isotope","_",batch))
            peptide <- exprs(peptide_batch) |> subset(rownames(exprs(peptide_batch)) %in% isotope)
            peptide <- t(peptide)
            colnames(peptide) <- paste0("X",colnames(peptide))
            vep <- read.delim(paste0("~/Caprion/analysis/peptide/",p,"/","vep.txt")) %>%
                   setNames(c("isotope_SNP", names(.)[-1])) %>%
                   extract(col = 1, into = c("isotope", "SNP"), regex = "^([0-9]+).tab:(.*)$") %>%
                   mutate(isotope=as.numeric(isotope)) %>%
                   select(isotope,SNP,Feature,Feature_type,Consequence,Codons,IMPACT,DISTANCE,SYMBOL,NEAREST)
            signals <- read.csv(paste0("~/Caprion/analysis/peptide/",p,"/",p,".cis.vs.trans"),header=TRUE) %>%
                       select(-geneStart,-geneEnd,-SNPChrom,-SNPPos) %>%
                       left_join(vep)
            for (type in c("protein","dr"))
            {
               prot <- get(type)
               col_sel <- !is.na(apply(peptide,2,sum))
               peptide_x <- peptide[,col_sel]
               fit <- glmnet(peptide_x,prot)
               plot(fit,label=TRUE)
               title(paste(paste0(p,"_",type),batch,"elastic net"),line=2.5)
               coef(fit, s = 0.1)
               beta <- coef(fit, s = 0.1)
               selected <- (beta!=0)
               selected_groups <- gsub("^X", "",rownames(beta)[selected[-1]])
               cvfit <- cv.glmnet(peptide_x, prot)
               plot(cvfit)
               title(paste(paste0(p,"_",type),batch,"CV elastic net"),line=2.5)
               coef(cvfit, s = "lambda.min")
               if (length(selected_groups)>0)
               {
                  cat("\n",p,type,"Peptide contributions:",selected_groups,"\n")
                  selected_signals <- subset(signals,isotope%in%selected_groups)
                  print(knitr::kable(selected_signals,caption=paste(p,batch)))
                  cat("\nCV min(lamda) =",cvfit$lambda.min,"\n")
               }
               if (TRUE)
               {
                 par(mfrow=c(2,2))
                 set.seed(12345)
                 fit <- adapense_cv(peptide,prot,alpha=0.75,cv_k=5,cv_repl=50,cl=cluster)
                 summary(fit,lambda="se")
                 set.seed(12345)
                 fit_075_1 <- adapense_cv(peptide,prot,alpha=0.75,exponent=1,cv_k=5,cv_repl=50,cl=cluster)
                 plot(fit_075_1)
                 set.seed(12345)
                 fit_100_1 <- adapense_cv(peptide,prot,alpha=1.00,exponent=1,cv_k=5,cv_repl=50,cl=cluster)
                 plot(fit_100_1)
                 set.seed(12345)
                 fit_075_2 <- adapense_cv(peptide,prot,alpha=0.75,exponent=2,cv_k=5,cv_repl=50,cl=cluster)
                 plot(fit_075_2)
                 set.seed(12345)
                 fit_100_2 <- adapense_cv(peptide,prot,alpha=1.00,exponent=2,cv_k=5,cv_repl=50,cl=cluster)
                 plot(fit_100_2)
                 prediction_performance(fit_075_1,fit_100_1,fit_075_2,fit_100_2,lambda="min")
               }
               par(mfrow=c(2,1))
            }
          }, error = function(e) {
            print(paste("Error at index", batch, ":", e$message))
          })
       }
       sink()
       par(opar)
       dev.off()
    }
  '
}

function glmnet_vignette()
{
  Rscript -e '
    library(glmnet)
  # glmnet
    data(QuickStartExample)
    x <- QuickStartExample$x
    y <- QuickStartExample$y
    fit <- glmnet(x, y)
    plot(fit,label=TRUE)
    print(fit)
    coef(fit, s = 0.1)
  # cv.glmnet
    set.seed(29)
    nx <- matrix(rnorm(5 * 20), 5, 20)
    predict(fit, newx = nx, s = c(0.1, 0.05))
    cvfit <- cv.glmnet(x, y)
    plot(cvfit)
    cvfit$lambda.min
    coef(cvfit, s = "lambda.min")
    predict(cvfit, newx = x[1:5,], s = "lambda.min")
  # weights
    wts <-  c(rep(1,50), rep(2,50))
    fit <- glmnet(x, y, alpha = 0.2, weights = wts, nlambda = 20)
    print(fit)
    fit <- glmnet(x, y)
    coef.apprx <- coef(fit, s = 0.5, exact = FALSE)
    coef.exact <- coef(fit, s = 0.5, exact = TRUE, x=x, y=y)
    cbind2(coef.exact[which(coef.exact != 0)],
           coef.apprx[which(coef.apprx != 0)])
    predict(fit, newx = x[1:5,], type = "response", s = 0.05)
    plot(fit, xvar = "lambda", label = TRUE)
    plot(fit, xvar = "dev", label = TRUE)
    cvfit <- cv.glmnet(x, y, type.measure = "mse", nfolds = 20)
    print(cvfit)
    structure(c(2.44, 0.08, 2.518, 0, 0), class = "proc_time", .Names = c("user.self",
    "sys.self", "elapsed", "user.child", "sys.child"))
    structure(c(0.508999999999999, 0.057, 1.56699999999999, 1.941,
    0.1), class = "proc_time", .Names = c("user.self", "sys.self",
    "elapsed", "user.child", "sys.child"))
    cvfit$lambda.min
    predict(cvfit, newx = x[1:5,], s = "lambda.min")
    coef(cvfit, s = "lambda.min")
    foldid <- sample(1:10, size = length(y), replace = TRUE)
    cv1  <- cv.glmnet(x, y, foldid = foldid, alpha = 1)
    cv.5 <- cv.glmnet(x, y, foldid = foldid, alpha = 0.5)
    cv0  <- cv.glmnet(x, y, foldid = foldid, alpha = 0)
    par(mfrow = c(2,2))
    plot(cv1); plot(cv.5); plot(cv0)
    plot(log(cv1$lambda)   , cv1$cvm , pch = 19, col = "red",
         xlab = "log(Lambda)", ylab = cv1$name)
    points(log(cv.5$lambda), cv.5$cvm, pch = 19, col = "grey")
    points(log(cv0$lambda) , cv0$cvm , pch = 19, col = "blue")
    legend("topleft", legend = c("alpha= 1", "alpha= .5", "alpha 0"),
           pch = 19, col = c("red","grey","blue"))
    tfit <- glmnet(x, y, lower.limits = -0.7, upper.limits = 0.5)
    plot(tfit)
    p.fac <- rep(1, 20)
    p.fac[c(1, 3, 5)] <- 0
    pfit <- glmnet(x, y, penalty.factor = p.fac)
    plot(pfit, label = TRUE)
  # parallel
    library(doMC)
    registerDoMC(cores = 2)
    X <- matrix(rnorm(1e4 * 200), 1e4, 200)
    Y <- rnorm(1e4)
    system.time(cv.glmnet(X, Y))
    system.time(cv.glmnet(X, Y, parallel = TRUE))
  '
}

glmnet_pense
