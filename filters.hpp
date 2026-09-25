#ifndef FILTERS_HPP
#define FILTERS_HPP

#include "dsptypes.hpp"

namespace dsp {

/* FIR filter: y[n] = sum_{k=0}^{num_taps-1} h[k] * x[n-k], for n in [0, sz).
 * Samples before x[0] are treated as absent (not zero-padded); taps beyond
 * the available history are simply skipped. x and y may be the same length
 * (sz); h has num_taps elements. */
void fir_filter(const Complex16* x, const Complex16* h, Complex16* y, size_t sz, size_t num_taps);

/* Direct-form biquad: y[t] = b0*x[t] + b1*x[t-1] + b2*x[t-2] - a1*y[t-1] - a2*y[t-2],
 * for t in [2, sz). y[0] and y[1] are left untouched by this function and must
 * be initialized by the caller. */
void biquad_filter(const Complex16* x, Q15 b0, Q15 b1, Q15 b2, Q15 a1, Q15 a2, Complex16* y, size_t sz);

/* Downsample by M: y[i] = x[M*i], for i in [0, sz/M). y must hold sz/M
 * elements. No anti-alias filtering is applied; filter first if needed. */
void decimate(const Complex16* x, Complex16* y, size_t M, size_t sz);

/* Upsample by M via zero-stuffing: for each input sample, M-1 zeros
 * followed by the sample itself. x has sz elements, y must hold sz*M.
 * No interpolation filtering (image rejection) is applied; filter the
 * result if needed. */
void interpolate(const Complex16* x, Complex16* y, size_t M, size_t sz);

} // namespace dsp

#endif // FILTERS_HPP
