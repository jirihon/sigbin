/**
 * Conversion of DNA sequence into numeric signal.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/11/03
 * Package: sigbin
 */

#ifndef SIGBIN_SIGNAL_HEADER
#define SIGBIN_SIGNAL_HEADER

/**
 * Convert sequence into numeric signal.
 *
 * @param seq DNA sequence.
 * @param signal Reference to output signal vector.
 */
template <typename T>
void sequence_to_signal(const std::string &seq, T signal) {
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

#endif
