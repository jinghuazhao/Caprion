## PCA-autoencoder comparison

The example compares performances of PCA and autoencoder for dimension reduction, which is adapted from 
https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/ and shown 
in[pca_ae_test.html](pca_ae_test.html), so that `the autoencoder is better at reconstructing the original data set than PCA 
when k is small, however the error converges as k increases. For very large data sets this difference will be larger and 
means a smaller data set could be used for the same error as PCA. When dealing with big data this is an important property`.

**REFERENCES (Variatinoal AutoEncoder)**

Hofert M, Prasad A, Zhu M (2019). Quasi-Monte Carlo for multivariate distributions viagenerative neural networks. https://arxiv.org/abs/1811.00683, https://CRAN.R-project.org/package=gnn.

Kingma DP, Welling M (2014). Auto-Encoding Variational Bayes. https://arxiv.org/abs/1312.6114, https://keras.rstudio.com/articles/examples/variational_autoencoder.html.

Trivadis SK (2017). Variational autoencoders for anomaly detection. https://rpubs.com/zkajdan/308801.
