###
## Signal-based binning.
##
## Author: Jiri Hon <jiri.hon@gmail.com>
## Date: 2016/11/03
## Package: sigbin
##

.sig_base_values = c(
  A = pi/4,
  C = -3*pi/4,
  G = 3*pi/4,
  T = -pi/4,
  N = 0
)

#' Compute Hjorth parameters for DNA sequences or numeric signal.
#'
#' @param x DNAString, DNAStringSet, numeric vector or list of numeric vectors.
#' @return Data frame of Hjorth parameters.
#'
#' @examples
#' sig_list <- dna_signal(DNAStringSet(c("ACGT","ACCT","GAAC")))
#' hp <- hjorth_params(sig_list)
#'
hjorth_params <- function(x) {
  if (class(x) == "numeric") {
    params <- vector("numeric", 3)
    hjorth_params_single_sig_cpp(x, params)
    return(data.frame(
      activity = params[1],
      mobility = params[2],
      complexity = params[3]))
  } else if (is.list(x)) {
    activity <- vector("numeric", length(x))
    mobility <- vector("numeric", length(x))
    complexity <- vector("numeric", length(x))
    hjorth_params_multi_sig_cpp(x, activity, mobility, complexity)
    return(data.frame(
      activity = activity,
      mobility = mobility,
      complexity = complexity))
  } else {
    stop("Signal must be eighter a numeric vector or a list of numeric vectors.")
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
    print(read_status)
  }
  closeAllConnections()
  return(data.frame(
    activity = activity,
    mobility = mobility,
    complexity = complexity)
  )
}

#' Convert DNA sequence into numeric signal.
#'
#' @param subject DNAString or DNAStringSet object.
#' @return Numeric vector with signal.
#'
#' @examples
#' sig <- dna2signal(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' sig_list <- dna2signal(DNAStringSet(c("ACGT","ACCT","GAAC")))
#'
dna_signal <- function(subject) {
  if (class(subject) == "DNAString") {
    return(sequence_to_signal_cpp(as.character(subject)))
  } else if (class(subject) == "DNAStringSet") {
    signal_list <- vector("list", length(subject))
    seq_set <- as.character(subject)
    for (i in 1:length(seq_set)) {
      signal_list[[i]] <- sequence_to_signal_cpp(seq_set[i])
    }
    return(signal_list)
  } else {
    stop("Sequence must be either DNAString or DNAStringSet object.")
  }
  return(NULL)
}
