/**
 * Computation of Hjorth parameters
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
 * Variance with Bessel's correction.
 *
 * @param sum_x Sum of all elements.
 * @param sum_x2 Sum of all squared elements.
 * @param n Number of elements
 */
inline double bvar(double sum_x, double sum_x2, double n) {
  double mean_x = sum_x / n;
  return (sum_x2 - 2*mean_x*sum_x + mean_x*mean_x*n) / (n - 1);
}

/**
 * Template for computing Hjorth parameters from generic signal.
 *
 * @param signal Numeric signal.
 * @param len Signal length.
 * @param activity Activity output array.
 * @param mobility Mobility output array.
 * @param complexity Complexity output array.
 */
template <typename T> void hjorth_params(T signal, int len, double *activity, double *mobility, double *complexity) {
  double prev_x = 0, prev_dx = 0,
    sum_x = 0, sum_x2 = 0,
    sum_dx = 0, sum_dx2 = 0,
    sum_ddx = 0, sum_ddx2 = 0;

  for (int i = 0; i < len; ++i) {
    double dx = signal[i] - prev_x;
    double ddx = dx - prev_dx;

    sum_x += signal[i];
    sum_x2 += signal[i]*signal[i];
    sum_dx += dx;
    sum_dx2 += dx*dx;
    sum_ddx += ddx;
    sum_ddx2 += ddx*ddx;

    prev_x = signal[i];
    prev_dx = dx;
  }
  *activity = bvar(sum_x, sum_x2, len);
  double sd_x = sqrt(*activity);
  double sd_dx = sqrt(bvar(sum_dx, sum_dx2, len));
  *mobility = sd_dx / sd_x;
  *complexity = (sqrt(bvar(sum_ddx, sum_ddx2, len)) / sd_dx) / *mobility;
}

/**
 * Compute Hjorth parameters for all sequences from character vector.
 *
 * @param seq_set Character vector containing set of sequences.
 * @param activity Activity output vector.
 * @param mobility Mobility output vector.
 * @param complexity Complexity output vector.
 * @param offset Index offset for writing to output vectors.
 */
// [[Rcpp::export]]
void hjorth_params_seq_cpp(
    std::vector<std::string> &seq_set,
    NumericVector &activity,
    NumericVector &mobility,
    NumericVector &complexity,
    int offset)
{
  int max_seq_len = 0;
  for (vector<string>::iterator it = seq_set.begin(); it != seq_set.end(); ++it) {
    if ((*it).length() > max_seq_len)
      max_seq_len = (*it).length();
  }
  double *signal = new double[max_seq_len];
  double act, mob, com;

  for (vector<string>::iterator it = seq_set.begin(); it != seq_set.end(); ++it) {
    string &seq = *it;
    sequence_to_signal<double *>(seq, signal);
    hjorth_params<double *>(signal, seq.length(), &act, &mob, &com);
    activity[offset] = act;
    mobility[offset] = mob;
    complexity[offset] = com;
    ++offset;
  }
  delete signal;
}

// [[Rcpp::export]]
void hjorth_params_single_sig_cpp(NumericVector &signal, NumericVector &result)
{
  double act, mob, com;
  hjorth_params<NumericVector &>(signal, signal.length(), &act, &mob, &com);
  result[0] = act;
  result[1] = mob;
  result[2] = com;
}

/**
 * Compute Hjorth parameters for all signals from a list.
 *
 * @param sig_list List containing set of numeric vectors.
 * @param activity Activity output vector.
 * @param mobility Mobility output vector.
 * @param complexity Complexity output vector.
 */
// [[Rcpp::export]]
void hjorth_params_sig_cpp(
    List &sig_list,
    NumericVector &activity,
    NumericVector &mobility,
    NumericVector &complexity)
{
  double act, mob, com;
  for (int i = 0; i < sig_list.length(); ++i) {
    NumericVector signal(sig_list[i]);
    hjorth_params<NumericVector &>(signal, signal.length(), &act, &mob, &com);
    activity[i] = act;
    mobility[i] = mob;
    complexity[i] = com;
  }
}
