#ifndef FFT_HPP
#define FFT_HPP

#include "dsptypes.hpp"

namespace dsp {

/* Radix-2 Cooley-Tukey FFT (decimation in time). N must be a power of two,
 * and at most 0x10000 (the phase resolution of the CORDIC twiddle factors).
 * x and X each hold N elements; they may safely alias the same buffer for
 * an in-place transform (x is fully consumed before X is written).
 *
 * inv:  false = forward transform (twiddle e^{-j2*pi*k/N}),
 *       true  = inverse transform (twiddle e^{+j2*pi*k/N}).
 *
 * Twiddle factors are computed once per call via sin_cos_q15() (see
 * dspmath.hpp for set_sin_cos_mode() to pick CORDIC vs. LUT) and reused
 * across every recursion level, not recomputed per stage.
 *
 * scale_mode:
 *   's' = scale every stage's output by 1/2, for an overall 1/N scale on
 *         the result. Keeps intermediate values within Q15 range
 *         regardless of input amplitude or direction. (default)
 *   'u' = unscaled, full-precision output. Not yet implemented.
 */
void fft(Complex16* x, size_t N, Complex16* X, char scale_mode = 's', bool inv = false);

} // namespace dsp

#endif // FFT_HPP
