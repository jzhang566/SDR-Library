#ifndef ANALOG_HPP
#define ANALOG_HPP

#include "dsptypes.hpp"

namespace dsp {

/* AM and SSB, as an NCO-driven mix onto an explicit carrier (a phase_t
 * step per sample), not zero-IF baseband -- SSB in particular has no
 * meaningful "baseband" form (it's defined by placing a Hilbert-shifted
 * signal onto a carrier), and AM follows the same shape for consistency.
 * start_phase/carrier_step let a caller continue a running NCO phase
 * across multiple calls if generating a long waveform in chunks (this
 * library's usual burst/block model: no state is kept internally). */

/* out[i] = (0.5 + 0.5*mod_index*audio[i]) * carrier[i]. The 0.5 bias/
 * scale is a Q15 range-safety measure -- the physical AM envelope
 * (1 + mod_index*audio) can reach 2.0 at full modulation, which doesn't
 * fit Q15's [-1, ~1) range at all. am_demodulate() recovers this same
 * 0.5-scaled, DC-biased envelope; remove the DC (e.g. via biquad_filter
 * as a simple high-pass) and rescale by 2/mod_index for calibrated audio
 * back, if needed. */
void am_modulate(const Q15* audio, size_t n, Q15 mod_index, phase_t start_phase, phase_t carrier_step, Complex16* out);

/* Envelope detection: sqrt_q15(x[i].magsq()). */
void am_demodulate(const Complex16* x, size_t n, Q15* audio_out);

/* SSB, phasing method. Unlike a classic phasing-method implementation,
 * which produces a real bandpass RF waveform, this produces a complex
 * single-sideband signal directly: out[i] = (audio +/- j*hilbert(audio))
 * * carrier[i], consistent with the rest of this library's complex-
 * baseband representation (AM, QPSK/QAM, ...). Taking .re of the output
 * recovers the classic real phasing-method formula exactly, if a real
 * waveform is needed downstream.
 *
 * audio and its Hilbert transform are combined into a complex "analytic"
 * signal before the carrier multiply, which has the same vector-
 * magnitude-vs-per-axis-amplitude issue QPSK/QAM's constellation points
 * have (see modulation.cpp's axis_step()): the combined vector magnitude
 * can be up to sqrt(2)x either input's amplitude, and a carrier rotation
 * can land that full magnitude on one axis. This scales the analytic
 * signal down by the same conservative ~1/sqrt(2) approximation used
 * there, so audio and its Hilbert transform should individually use at
 * most Q15's full range -- the scale-down is handled internally, not
 * something the caller needs to pre-apply. ssb_demodulate() recovers
 * this same scaled amplitude.
 *
 * hilbert_taps must be odd (see design_hilbert_q15()). */
void ssb_modulate(const Q15* audio, size_t n, size_t hilbert_taps, bool upper_sideband,
                   phase_t start_phase, phase_t carrier_step, Complex16* out);

/* Mixes down by the carrier and takes the real part -- for this
 * library's own complex SSB representation (see ssb_modulate()), no
 * further Hilbert transform is needed on receive. */
void ssb_demodulate(const Complex16* x, size_t n, phase_t start_phase, phase_t carrier_step, Q15* audio_out);

} // namespace dsp

#endif // ANALOG_HPP
