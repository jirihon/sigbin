###
## Hjorth parameters.
##
## Author: Jiri Hon <jiri.hon@gmail.com>
## Date: 2016/11/03
## Package: sigbin
##

# Compute Hjorth parameters for single signal vector.
#
# @param x Signal vector.
# @return Data frame of Hjorth parameters.
#
.hjorth_params_single_sig <- function(x) {
  params <- vector("numeric", 3)
  hjorth_params_single_sig_cpp(x, params)
  return(data.frame(
    activity = params[1],
    mobility = params[2],
    complexity = params[3])
  )
}

# Compute Hjorth parameters for multiple signals.
#
# @param x Signal vector.
# @param fn Computation function.
# @param args Additional arguments for computation function.
# @return Data fram of Hjorth parameters.
#
.hjorth_params_multi <- function(x, fn, args = list()) {
  activity <- vector("numeric", length(x))
  mobility <- vector("numeric", length(x))
  complexity <- vector("numeric", length(x))
  do.call(fn, unlist(list(list(x), args, list(activity, mobility, complexity)), recursive = FALSE))
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
    return(.hjorth_params_single_sig(x))
  } else if (is.list(x)) {
    return(.hjorth_params_multi(x, hjorth_params_multi_sig_cpp))
  } else if (class(x) == "DNAString") {
    return(.hjorth_params_single_sig(sequence_to_signal_cpp(as.character(x))))
  } else if (class(x) == "DNAStringSet") {
    seq_set <- as.character(x)
    max_seq_len <- max(sapply(seq_set, nchar))
    return(.hjorth_params_multi(seq_set, hjorth_params_multi_seq_cpp, list(max_seq_len, 0)))
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
#'
fasta_hjorth_params <- function(filepath, block_size = 1e4) {
  seqlens <- fasta.seqlengths(filepath)
  activity <- vector("numeric", length(seqlens))
  mobility <- vector("numeric", length(seqlens))
  complexity <- vector("numeric", length(seqlens))
  max_seq_len <- max(seqlens)

  read_status <- 0
  while (read_status < length(seqlens)) {
    dnaset <- readDNAStringSet(filepath, nrec = block_size, skip = read_status, use.names = FALSE)
    seqset <- as.character(dnaset)
    hjorth_params_multi_seq_cpp(seqset, max_seq_len, read_status, activity, mobility, complexity)
    read_status <- read_status + length(dnaset)

    status <- ceiling(read_status / length(seqlens) * 100)
    cat(sprintf("\rStatus: %d %%", status))
    flush(stdout())
  }
  return(data.frame(
    activity = activity,
    mobility = mobility,
    complexity = complexity)
  )
}
