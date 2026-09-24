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

} // namespace dsp

#endif // DSP_MATH_HPP
