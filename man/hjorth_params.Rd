% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hjorth.R
\name{hjorth_params}
\alias{hjorth_params}
\title{Compute Hjorth parameters for DNA sequences or numeric signal.}
\usage{
hjorth_params(x)
}
\arguments{
\item{x}{DNAString, DNAStringSet, numeric vector or list of numeric vectors.}
}
\value{
Data frame of Hjorth parameters.
}
\description{
Compute Hjorth parameters for DNA sequences or numeric signal.
}
\examples{
dna_set <- DNAStringSet(c("ACGT","ACCT","GAAC"))
hp1 <- hjorth_params(dna_set)
hp2 <- hjorth_params(dna_set[[1]])

sig_list <- dna_signal(dna_set)
hp3 <- hjorth_params(sig_list)
hp4 <- hjorth_params(sig_list[[1]])

}

