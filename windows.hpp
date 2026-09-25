#ifndef WINDOWS_HPP
#define WINDOWS_HPP

#include "dsptypes.hpp"

namespace dsp {

/* All windows below apply an sz-point window to the first sz samples of x
 * in place, then zero-pad samples [sz, N) -- so x ends up ready for an
 * N-point FFT in one call. sz must be in (0, N]; if sz == N (the common
 * case) there's nothing to zero-pad. */

/* Hamming window: w[n] = 0.54 - 0.46*cos(2*pi*n/(sz-1)). */
void hamming_window(Complex16* x, size_t sz, size_t N);

/* 4-term minimum Blackman-Harris window (~-92 dB sidelobes):
 * w[n] = 0.35875 - 0.48829*cos(t) + 0.14128*cos(2t) - 0.01168*cos(3t),
 * t = 2*pi*n/(sz-1). */
void blackman_harris_window(Complex16* x, size_t sz, size_t N);

/* Kaiser window with shape parameter beta (larger beta = wider mainlobe,
 * lower sidelobes; beta ~= 8.6 gives roughly -100 dB sidelobes):
 * w[n] = I0(beta*sqrt(1 - r^2)) / I0(beta), r = 2n/(sz-1) - 1, where I0 is
 * the zeroth-order modified Bessel function of the first kind. I0(beta)
 * can be far outside Q15 range for typical beta, so unlike the other
 * windows this one computes its Bessel series in double (a one-time
 * setup cost, not a per-sample operation) and quantizes the final
 * (bounded, [0,1]) ratio down to Q15. */
void kaiser_window(Complex16* x, size_t sz, size_t N, double beta);

} // namespace dsp

#endif // WINDOWS_HPP
