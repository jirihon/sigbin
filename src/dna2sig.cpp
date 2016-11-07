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


template <typename T> void sequence_to_signal(const string &seq, T signal) {
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
}

inline double bvar(double sum_x, double sum_x2, double n) {
  double mean_x = sum_x / n;
  return (sum_x2 - 2*mean_x*sum_x + mean_x*mean_x*n) / (n - 1);
}

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

// [[Rcpp::export]]
void hjorth_params_multi_seq_cpp(
    std::vector<std::string> &seq_set,
    int max_seq_len,
    int offset,
    NumericVector &activity,
    NumericVector &mobility,
    NumericVector &complexity)
{
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

// [[Rcpp::export]]
void hjorth_params_multi_sig_cpp(
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
