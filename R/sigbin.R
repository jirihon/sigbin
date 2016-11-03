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
  T = -pi/4)


.sequence_signal <- function(dnastring) {
  return(toComplex(dnastring, .sig_base_values))
}

.activity <- function(x) {
  var(x)
}

.mobility <- function(x, dx) {
  sd(dx)/sd(x)
}

.complexity <- function(x, dx) {
  ddx <- diff(c(0, dx))
  (sd(ddx)/sd(dx)) / (sd(dx)/sd(x))
}

#' Compute Hjorth parameters from numeric signal.
#'
#' @param x Numeric vector or list of numeric vectors.
#' @return Data frame of Hjorth parameters.
#'
#' @examples
#' sig_list <- dna2signal(DNAStringSet(c("ACGT","ACCT","GAAC")))
#' hp <- hjorth_params(sig_list)
#'
hjorth_params <- function(x) {
  if (class(x) == "numeric") {
    dx <- diff(c(0, x))
    return(data.frame(
      activity = .activity(x),
      mobility = .mobility(x, dx),
      complexity = .complexity(x, dx)))
  } else if (class(x) == "list") {
    activity <- vector("numeric", length(x))
    mobility <- vector("numeric", length(x))
    complexity <- vector("numeric", length(x))

    for (i in 1:length(x)) {
      if (class(x[[i]]) == "numeric") {
        dx <- diff(c(0, x[[i]]))
        activity[i] <- .activity(x[[i]])
        mobility[i] <- .mobility(x[[i]], dx)
        complexity[i] <- .complexity(x[[i]], dx)
      } else {
        stop("All elements of signal list must be numeric vectors.")
      }
    }
    return(data.frame(
      activity = activity,
      mobility = mobility,
      complexity = complexity))
  } else {
    stop("Signal must be eighter numeric vector or list of numeric vectors.")
  }
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
    return(.sequence_signal(subject))
  } else if (class(subject) == "DNAStringSet") {
    signal_list <- vector("list", length(subject))
    for (i in 1:length(subject)) {
      signal_list[[i]] <- .sequence_signal(subject[[i]])
    }
    return(signal_list)
  } else {
    stop("Sequence must be either DNAString or DNAStringSet object.")
  }
  return(NULL)
}
