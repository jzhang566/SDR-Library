#ifndef DSP_MATH_HPP
#define DSP_MATH_HPP

#include "dsptypes.hpp"

namespace dsp {

/* ===== Phase / complex helpers ===== */

/* Phase (argument) of a complex value via CORDIC. */
phase_t phase_c16(const Complex16& a);

/* ===== CORDIC ===== */

/* Compute sine and cosine of a phase using CORDIC, as the unit phasor
 * (re = cos, im = sin). */
Complex16 sin_cos_cordic_q15(phase_t phase);

/* Compute the phase (arctangent) of a complex vector using CORDIC. */
phase_t atan_cordic_q15(const Complex16& v);

/* ===== Sinusoid lookup table =====
 *
 * ptr is a user-allocated buffer of size (sz * sizeof(Q15)).
 * sz must be a power of 2, greater than 4.
 * Only one sin-cos LUT can be created at a time.
 * Returns true on success, false on invalid size.
 */
bool generate_sinusoid_lut_q15(Q15* ptr, size_t sz);

/* Cosine via LUT with linear interpolation. */
Q15 cos_lut_q15(phase_t angle);

/* Sine via LUT with linear interpolation. */
Q15 sin_lut_q15(phase_t angle);

/* Cosine via LUT, nearest-neighbour (no interpolation). */
Q15 cos_lut_q15_nointerp(phase_t angle);

/* Sine via LUT, nearest-neighbour (no interpolation). */
Q15 sin_lut_q15_nointerp(phase_t angle);

/* ===== Buffer / vector operations =====
 *
 * All take element count n and operate on out[i] = f(a[i], b[i]) for
 * i in [0, n). out may alias a and/or b.
 */

/* out[i] = a[i] + b[i] */
void pointwise_add(const Complex16* a, const Complex16* b, Complex16* out, size_t n);

/* out[i] = a[i] - b[i] */
void pointwise_sub(const Complex16* a, const Complex16* b, Complex16* out, size_t n);

/* out[i] = a[i] * b[i] */
void pointwise_mul(const Complex16* a, const Complex16* b, Complex16* out, size_t n);

/* out[i] = a[i] * b[i], b real-valued (e.g. applying a window function). */
void pointwise_mul(const Complex16* a, const Q15* b, Complex16* out, size_t n);

/* out[i] = a[i] * factor (e.g. gain/normalization). */
void scale(const Complex16* a, Q15 factor, Complex16* out, size_t n);

/* out[i] = |a[i]|^2 */
void magsq(const Complex16* a, Q15* out, size_t n);

/* Complex (Hermitian) inner product: sum_i a[i] * conj(b[i]).
 * Accumulates internally in 64-bit fixed point so it doesn't saturate
 * per-term on long vectors; the final sum is saturated to Q15 once. */
Complex16 dot(const Complex16* a, const Complex16* b, size_t n);

/* Signal energy: sum_i |a[i]|^2, equivalent to dot(a, a).re but computed
 * directly. Same wide accumulation as dot(); still saturates to Q15 on
 * return, so it's most useful over short buffers or pre-scaled inputs. */
Q15 energy(const Complex16* a, size_t n);

} // namespace dsp

#endif // DSP_MATH_HPP
