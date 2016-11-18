###
## Hjorth parameters.
##
## Author: Jiri Hon <jiri.hon@gmail.com>
## Date: 2016/11/03
## Package: sigbin
##

# Compute Hjorth parameters for a set of signals/sequences.
#
# @param x Set of signals/sequence.
# @param fn Backend computation function.
# @param args Additional arguments for computation function.
# @return Data frame of Hjorth parameters.
#
.hjorth_params <- function(x, fn, args = list()) {
  activity <- vector("numeric", length(x))
  mobility <- vector("numeric", length(x))
  complexity <- vector("numeric", length(x))
  do.call(fn, unlist(list(list(x, activity, mobility, complexity), args), recursive = FALSE))
  return(data.frame(
    activity = activity,
    mobility = mobility,
    complexity = complexity)
  )
}

#' Compute Hjorth parameters for DNA sequences or numeric signal.
#'
#' @param x DNAString, DNAStringSet, numeric vector or list of numeric vectors.
#' @return Data frame of Hjorth parameters.
#'
#' @examples
#' dna_set <- DNAStringSet(c("ACGT","ACCT","GAAC"))
#' hp1 <- hjorth_params(dna_set)
#' hp2 <- hjorth_params(dna_set[[1]])
#'
#' sig_list <- dna_signal(dna_set)
#' hp3 <- hjorth_params(sig_list)
#' hp4 <- hjorth_params(sig_list[[1]])
#'
hjorth_params <- function(x) {
  if (is.numeric(x)) {
    return(.hjorth_params(list(x), hjorth_params_sig_cpp))
  } else if (is.list(x)) {
    return(.hjorth_params(x, hjorth_params_sig_cpp))
  } else if (class(x) == "DNAString") {
    return(.hjorth_params(list(sequence_to_signal_cpp(as.character(x))), hjorth_params_sig_cpp))
  } else if (class(x) == "DNAStringSet") {
    seq_set <- as.character(x)
    max_seq_len <- max(sapply(seq_set, nchar))
    return(.hjorth_params(seq_set, hjorth_params_seq_cpp, list(0)))
  } else {
    stop("Signal must be eighter a numeric vector, list of numeric vectors, DNAString or DNAStringSet object.")
  }
}

#' Compute Hjorth parameters for all DNA sequences from FASTA file.
#'
#' @param filepath Path to FASTA file.
#' @param block_size Number of sequences read at once into memory.
#' @return Data frame of Hjorth parameters.
#'
#' @examples
#' file <- system.file("extdata", "reads.fasta", package="sigbin")
#' hp <- fasta_hjorth_params(file)
#'
fasta_hjorth_params <- function(filepath, block_size = 1e4) {
  seq_lens <- fasta.seqlengths(filepath)
  activity <- vector("numeric", length(seq_lens))
  mobility <- vector("numeric", length(seq_lens))
  complexity <- vector("numeric", length(seq_lens))

  read_status <- 0
  while (read_status < length(seq_lens)) {
    dna_set <- readDNAStringSet(filepath, nrec = block_size, skip = read_status, use.names = FALSE)
    seq_set <- as.character(dna_set)
    hjorth_params_seq_cpp(seq_set, activity, mobility, complexity, read_status)
    read_status <- read_status + length(dna_set)

    status <- ceiling(read_status / length(seq_lens) * 100)
    cat(sprintf("\rStatus: %d %%", status))
    flush(stdout())
  }
  return(data.frame(
    activity = activity,
    mobility = mobility,
    complexity = complexity)
  )
}
