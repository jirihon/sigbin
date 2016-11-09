/**
 * Conversion of DNA sequence into numeric signal.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/11/03
 * Package: sigbin
 */

#include <Rcpp.h>
#include "signal.h"

using namespace Rcpp;
using namespace std;

/**
 * Convert DNA sequence into a numeric signal.
 *
 * @param seq DNA sequence string
 * @return Numeric signal
 */
// [[Rcpp::export]]
NumericVector sequence_to_signal_cpp(std::string &seq) {
  NumericVector signal(seq.size());
  sequence_to_signal<NumericVector &>(seq, signal);
  return signal;
}
