/**
 * Conversion of DNA sequence into numeric signal.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/11/03
 * Package: sigbin
 */

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


/**
 * Convert DNA sequence into a numeric signal.
 *
 * @param seq DNA sequence string
 * @return Numeric signal
 */
NumericVector sequence2signal(string seq) {
  NumericVector signal(seq.size());
  const double s_a = M_PI/4;
  const double s_c = -3*M_PI/4;
  const double s_g = 3*M_PI/4;
  const double s_t = -M_PI/4;

  for (int i = 0; i < seq.size(); ++i) {
    switch(seq[i]) {
    case 'A':
      signal[i] = s_a;
      break;
    case 'C':
      signal[i] = s_c;
      break;
    case 'G':
      signal[i] = s_g;
      break;
    case 'T':
      signal[i] = s_t;
      break;
    default:
      signal[i] = 0;
    break;
    }
  }
  return signal;
}

//' Convert DNA sequence into numeric signal.
//'
//' @param subject DNAString or DNAStringSet object.
//' @return Numeric vector with signal
//'
//' @examples
//' sig <- dna2signal(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
//' sig_list <- dna2signal(DNAStringSet(c("ACGT","ACCT","GAAC")))
//'
SEXP dna2signal(SEXP subject) {
  Function R_as_character("as.character");
  Function R_get_class("class");
  CharacterVector sequence_class = R_as_character(R_get_class(subject));
  SEXP result = R_NilValue;

  if (sequence_class[0] == "DNAString") {
    string seq = as<string>(R_as_character(subject));
    result = sequence2signal(seq);
  } else if (sequence_class[0] == "DNAStringSet") {
    List signal_list;
    Function R_subset("[[");
    Function R_length("length");

    int set_len = as<int>(R_length(subject));
    for (int i = 1; i <= set_len; ++i) {
      string seq = as<string>(R_as_character(R_subset(subject, i)));
      signal_list.push_back(sequence2signal(seq));
    }
    result = signal_list;
  } else {
    stop("Sequence must be either DNAString or DNAStringSet object.");
  }
  return result;
}

// Testing code examples
/*** R
library(Biostrings)
sig <- dna2signal(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
sig_list <- dna2signal(DNAStringSet(c("ACGT","ACCT","GAAC")))
*/
