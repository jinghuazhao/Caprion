# 20-11-2019 JHZ

source("caprion.inc")

# https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/

# autoencoder in keras
suppressPackageStartupMessages(library(keras))
library(ggplot2)
library(plotly)
library(DAAG)

# standardise
minmax <- function(x) (x - min(x))/(max(x) - min(x))
x_train <- apply(ais[,1:11], 2, minmax)
# x_train <- apply(pheno_protein[,-(1:8)], 2, minmax)
# PCA
pca <- prcomp(x_train)
# plot cumulative plot
qplot(x = 1:11, y = cumsum(pca$sdev)/sum(pca$sdev), geom = "line")
ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, col = ais$sex)) + geom_point()
pca_plotly <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color = ~ais$sex) %>% add_markers()
suppressPackageStartupMessages(library(keras))

# set training data
x_train <- as.matrix(x_train)

# set model
model <- keras_model_sequential()
model %>%
  layer_dense(units = 6, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
  layer_dense(units = 6, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

# view model layers
summary(model)
# compile model
model %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam"
)

# fit model
model %>% fit(
  x = x_train, 
  y = x_train, 
  epochs = 2000,
  verbose = 0
)

# evaluate the performance of the model
mse.ae2 <- evaluate(model, x_train, x_train)
mse.ae2

# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, x_train)

ggplot(data.frame(PC1 = intermediate_output[,1], PC2 = intermediate_output[,2]), aes(x = PC1, y = PC2, col = ais$sex)) + geom_point()
# PCA reconstruction
pca.recon <- function(pca, x, k){
  mu <- matrix(rep(pca$center, nrow(pca$x)), nrow = nrow(pca$x), byrow = T)
  recon <- pca$x[,1:k] %*% t(pca$rotation[,1:k]) + mu
  mse <- mean((recon - x)^2)
  return(list(x = recon, mse = mse))
}

xhat <- rep(NA, 10)
for(k in 1:10){
  xhat[k] <- pca.recon(pca, x_train, k)$mse
}

ae.mse <- rep(NA, 5)
for(k in 1:5){
  modelk <- keras_model_sequential()
  modelk %>%
    layer_dense(units = 6, activation = "tanh", input_shape = ncol(x_train)) %>%
    layer_dense(units = k, activation = "tanh", name = "bottleneck") %>%
    layer_dense(units = 6, activation = "tanh") %>%
    layer_dense(units = ncol(x_train))

  modelk %>% compile(
    loss = "mean_squared_error", 
    optimizer = "adam"
  )

  modelk %>% fit(
    x = x_train, 
    y = x_train, 
    epochs = 5000,
    verbose = 0
  )

  ae.mse[k] <- unname(evaluate(modelk, x_train, x_train))
}

df <- data.frame(k = c(1:10, 1:5), mse = c(xhat, ae.mse), method = c(rep("pca", 10), rep("autoencoder", 5)))
ggplot(df, aes(x = k, y = mse, col = method)) + geom_line()

# ae(pheno_protein[,-c(1:9)],hidden.layers=c(987,197,987))

# UMAP
library(uwot)
pdf("umap.pdf")
df_umap <- umap(df,n_neighbors = 5, learning_rate = 0.5, pca=50, spread = 3)
plot(df_umap,col="black",cex=0.4)
points(df_umap[(1:197)[idx11],],col="blue",cex=0.4)
points(df_umap[(1:197)[idx3],],col="red",cex=0.4)
title(paste("UMAP with n_neighbors=5, learning_rate=0.5, pca=50, spread=3"))
dev.off()

# PCA
ppc <- with(pheno_protein, prcomp(na.omit(pheno_protein[,-(1:8)]), rank=50, scale=TRUE))
write.table(data.frame(with(pheno_protein,caprion_id),with(ppc,x)),file="pca.txt",sep="\t",quote=FALSE,row.names=FALSE)
pdf("pca.pdf")
screeplot(ppc, npcs=20, type="lines", main="PCA screeplot")
idx11 <- with(pheno_protein,caprion_id)%in%eleven
idx3 <- with(pheno_protein,caprion_id)%in%three
with(ppc, {
  plot(x[,1:2], main="PC1 -- PC2", cex=0.5, pch=21)
  points(x[idx11,1],x[idx11,2], col="blue", cex=0.8, pch=10)
  points(x[idx3,1],x[idx3,2], col="red", cex=0.8, pch=16)
  legend("bottom", legend = c("lower detection", "lipemic"), box.lty = 0, cex = 0.8,
         col = c("red", "blue"), horiz = TRUE, inset = c(0, 1), xpd = TRUE, pch = c(10, 16))
})
biplot(ppc,cex=0.1)
title("biplot")
dev.off()

# RLE plot
pdf("rle.pdf", width=50, height=12)
par(mfrow=c(2,1))
makeRLEplot(d, log2.data=FALSE, groups=group, col.group=col.group, cex=0.3, showTitle=TRUE)
makeRLEplot(d, log2.data=FALSE, cex=0.3, showTitle=TRUE, title="Uncoloured relative log expression (RLE) plot")
dev.off()

# overlaps
protein_list <- within(protein_list,{Protein <- gsub("_HUMAN","",Protein)})
caprion_inf()

# correlation
with(pheno_protein,cor(crp,CRP,use="complete.obs"))
with(pheno_protein,cor(transf,TRFE,use="complete.obs"))

# cluster analysis
km <- kmeans(d,3)
with(km, table(cluster))
library(mclust)
mc <- Mclust(d)
summary(mc)
pdf("mc.pdf")
plot(mc, what = "BIC")
dev.off()

# regression
attach(pheno_protein)
options(width=500)
cat("protein","intercept","sex",sep="\t",file="sex.tsv")
cat("\n",append=TRUE,file="sex.tsv")
cat("protein","intercept","sex","age","bmi",sep="\t",file="lm.tsv")
cat("\n",append=TRUE,file="lm.tsv")
sapply(seq(1, ncol(df)), regfun)
detach(pheno_protein)
pdf("qq.pdf")
sex <- read.delim("sex.tsv",as.is=TRUE)
gap::qqunif(with(sex,sex),cex=0.4)
title("QQ plot for sex")
lm <- read.delim("lm.tsv",as.is=TRUE)
gap::qqunif(with(lm,bmi),cex=0.4)
title("QQ plot for BMI")
dev.off()

# EDA plots
df <- t1
rownames(df) <- gsub("ZWK","",rownames(df))
pdf("scatter-histogram-boxwhisker.pdf")
par(mfrow=c(3,1))
sapply(seq(1, ncol(df)), plotfun)
dev.off()

# Box-Whisker plot
pdf("box.pdf", width=50, height=12)
np <- ncol(df)
xtick <- seq(1, np, by=1)
boxplot(df,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(df), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("raw data")
sdf <- scale(df,scale=FALSE)
boxplot(sdf,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(sdf), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("mean-centred data")
ssdf <- scale(df,scale=TRUE)
boxplot(ssdf,xaxt="n",cex=0.2)
axis(side=1, at=xtick, tick=FALSE, labels=FALSE)
text(x=xtick, par("usr")[3], labels = colnames(ssdf), srt=90, pos=1, xpd = TRUE, cex=0.25)
title("mean-centred and scaled data")
dev.off()
