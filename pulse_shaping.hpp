#ifndef PULSE_SHAPING_HPP
#define PULSE_SHAPING_HPP

#include "dsptypes.hpp"

namespace dsp {

/* Designs a root-raised-cosine FIR filter, usable directly with
 * fir_filter(). num_taps should be odd (gives an exact integer group
 * delay of (num_taps-1)/2 samples, and a true center tap at t=0).
 * rolloff is in [0, 1] (0 = sinc/brick-wall, 1 = maximum excess
 * bandwidth). sps is samples per symbol (the oversampling factor).
 *
 * Generated via the closed-form RRC impulse response in double (a
 * one-time setup cost, not a per-sample operation -- same precedent as
 * kaiser_window()'s Bessel series), with the formula's two removable
 * singularities (at t=0 and t=+-T/(4*rolloff)) evaluated via their
 * analytic limits rather than a raw 0/0. Taps are peak-normalized (the
 * largest tap sits at the Q15 edge) so no rolloff/sps/num_taps
 * combination can saturate -- this is an amplitude-safety
 * normalization, not a specific unit-gain-at-symbol-samples
 * calibration; a caller chaining TX/RX RRC filters for exact power
 * accounting will want to handle that scaling itself. */
void design_rrc_q15(Q15* taps, size_t num_taps, double rolloff, size_t sps);

} // namespace dsp

#endif // PULSE_SHAPING_HPP
