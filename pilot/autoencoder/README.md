# Autoencoder

As shown at <a href="https://www.r-bloggers.com/pca-vs-autoencoders-for-dimensionality-reduction/">R-bloggers<img src="https://i0.wp.com/gradientdescending.com/wp-content/uploads/2018/07/reconstruction-1.png" width="560" height="400" align="right"></a>,

> autoencoder is better at reconstructing the original data set than PCA when k is small, 

where `k` corresponds to the number of principal components in PCA or bottleneck dimension in AE,

> however the error converges as k increases. For very large data sets this difference will be
> larger and means a smaller data set could be used for the same error as PCA. When dealing 
> with big data this is an important property`.

The local adoption is [ae_test.Rmd](ae_test.Rmd) which produces [ae_test.html](ae_test.html) and [ae_test.pdf](ae_test.pdf).

Additional work will be on variatinoal autoencoder (VAE) and denoising counterpart as indicated in the references below.

## REFERENCES

Bishop CM, Bishop H (2024), Deep learning: foundations and concepts, Springer International Publishing, [DOI: 10.1007/978-3-031-45468-4](https://link.springer.com/book/10.1007/978-3-031-45468-4).

Bludau I, Frank M, Dörig C. et al. Systematic detection of functional proteoform groups from bottom-up proteomic datasets. Nat Commun 12, 3810 (2021). [https://doi.org/10.1038/s41467-021-24030-x](https://doi.org/10.1038/s41467-021-24030-x).

Hofert M, Prasad A, Zhu M (2019). Quasi-Monte Carlo for multivariate distributions viagenerative neural networks. [https://arxiv.org/abs/1811.00683](https://arxiv.org/abs/1811.00683), [https://CRAN.R-project.org/package=gnn](https://CRAN.R-project.org/package=gnn).

Kingma DP, Welling M (2014). Auto-Encoding Variational Bayes. [https://arxiv.org/abs/1312.6114](https://arxiv.org/abs/1312.6114), [https://keras.rstudio.com/articles/examples/variational_autoencoder.html](https://keras.rstudio.com/articles/examples/variational_autoencoder.html).

Ng A. Sparse autoencoder, CS294A Lecture notes, <https://web.stanford.edu/class/archive/cs/cs294a/cs294a.1104/sparseAutoencoder.pdf>.

Sattarov T, Herurkar D, Hees J (2023). Explaining Anomalies using Denoising Autoencoders for Financial Tabular Data. Technical Report 2023-01. Deutsche Bundesban.

Trivadis SK (2017). Variational autoencoders for anomaly detection. [https://rpubs.com/zkajdan/308801](https://rpubs.com/zkajdan/308801).

## URLs

<https://github.com/diazale/gt-dimred>, <https://github.com/lmcinnes/umap> (<https://umap-learn.readthedocs.io/en/latest/>) and <https://keras.io/examples/timeseries/timeseries_anomaly_detection/>, <https://www.mathworks.com/help/deeplearning/ug/anomaly-detection-using-autoencoder-and-wavelets.html>, among others.
