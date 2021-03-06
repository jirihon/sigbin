% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hjorth.R
\name{fasta_hjorth_params}
\alias{fasta_hjorth_params}
\title{Compute Hjorth parameters for all DNA sequences from FASTA file.}
\usage{
fasta_hjorth_params(filepath, block_size = 10000)
}
\arguments{
\item{filepath}{Path to FASTA file.}

\item{block_size}{Number of sequences read at once into memory.}
}
\value{
Data frame of Hjorth parameters.
}
\description{
Compute Hjorth parameters for all DNA sequences from FASTA file.
}
\examples{
file <- system.file("extdata", "reads.fasta", package="sigbin")
hp <- fasta_hjorth_params(file)

}

