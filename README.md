# How to install

Download the latest stable R from <https://cloud.r-project.org/> (possibly optional).

Install Bioconductor and `sigbin` package dependencies:
```R
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "fastcluster", "dynamicTreeCut", "fields", "rgl", "Rcpp", "devtools"))
```

Install `sigbin` package directly from GitHub:
```R
devtools::install_github("jirihon/sigbin")
```

## Implemented functions

Following functions are implemented: `dna_signal`, `hjorth_params`, `fasta_hjorth_params`, `sigbin` and `plot_sigbin_3d`. Each function has a documentation entry and you can run its example code by `example` built-in function.
