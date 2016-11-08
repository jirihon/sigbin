###
## DNA to signal conversion.
##
## Author: Jiri Hon <jiri.hon@gmail.com>
## Date: 2016/11/03
## Package: sigbin
##

#' Convert DNA sequence into numeric signal.
#'
#' @param subject DNAString or DNAStringSet object.
#' @return Numeric vector with signal.
#'
#' @examples
#' sig <- dna_signal(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' sig_list <- dna_signal(DNAStringSet(c("ACGT","ACCT","GAAC")))
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
