# How to install

Download the latest stable R from <https://cloud.r-project.org/> (possibly optional).

Install Bioconductor and `sigbin` package dependencies:
```R
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "fastcluster", "dynamicTreeCut", "Rcpp", "devtools"))
```

Install `sigbin` package directly from GitHub:
```R
devtools::install_github("jirihon/sigbin")
```
