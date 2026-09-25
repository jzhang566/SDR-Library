#ifndef HILBERT_HPP
#define HILBERT_HPP

#include "dsptypes.hpp"

namespace dsp {

/* Designs a windowed FIR Hilbert transformer (broadband ~90-degree phase
 * shift), for the SSB phasing method. num_taps must be odd (Type III:
 * antisymmetric, zero center tap). Ideal (infinite) coefficients are
 * 2/(pi*m) for odd m, 0 for even m (m = tap index - center); windowed
 * with blackman_harris_window() to control Gibbs ringing.
 *
 * Unlike design_rrc_q15()/the window functions, this does NOT
 * peak-normalize to the Q15 edge: unity passband gain matters here (the
 * Hilbert-shifted signal must have the same amplitude as the original
 * for the phasing method's audio*cos -+ hilbert(audio)*sin combine to be
 * meaningful), so the natural ideal-response scale (peak ~0.6366,
 * comfortably within Q15 range) is kept, only tapered by the window. */
void design_hilbert_q15(Q15* taps, size_t num_taps);

} // namespace dsp

#endif // HILBERT_HPP
