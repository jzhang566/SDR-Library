#include "../analog.hpp"
#include "../dspmath.hpp"
#include "../filters.hpp"
#include "../hilbert.hpp"

#include <cassert>
#include <vector>

namespace dsp {

void am_modulate(const Q15* audio, size_t n, Q15 mod_index, phase_t start_phase, phase_t carrier_step, Complex16* out) {
    constexpr Q15 kHalf(static_cast<int16_t>(16384)); // 0.5
    phase_t phase = start_phase;
    for (size_t i = 0; i < n; ++i) {
        Q15 envelope = kHalf + kHalf * mod_index * audio[i];
        out[i] = sin_cos_q15(phase) * envelope;
        phase = static_cast<phase_t>(phase + carrier_step);
    }
}

void am_demodulate(const Complex16* x, size_t n, Q15* audio_out) {
    for (size_t i = 0; i < n; ++i) {
        audio_out[i] = sqrt_q15(x[i].magsq());
    }
}

void ssb_modulate(const Q15* audio, size_t n, size_t hilbert_taps, bool upper_sideband,
                   phase_t start_phase, phase_t carrier_step, Complex16* out) {
    assert(hilbert_taps % 2 == 1 && "ssb_modulate: hilbert_taps should be odd");

    std::vector<Q15> h(hilbert_taps);
    design_hilbert_q15(h.data(), hilbert_taps);
    std::vector<Complex16> h_c(hilbert_taps);
    for (size_t i = 0; i < hilbert_taps; ++i) h_c[i] = Complex16(h[i], Q15(0));

    std::vector<Complex16> audio_c(n);
    for (size_t i = 0; i < n; ++i) audio_c[i] = Complex16(audio[i], Q15(0));

    std::vector<Complex16> hilbert_out(n);
    fir_filter(audio_c.data(), h_c.data(), hilbert_out.data(), n, hilbert_taps);

    // Conservative 1/sqrt(2) approximation, same rationale as
    // modulation.cpp's axis_step(): bounds the combined (audio,
    // hilbert(audio)) vector magnitude so a subsequent carrier rotation
    // can't saturate either axis.
    constexpr Q15 kInvSqrt2(static_cast<int16_t>(22938)); // ~0.7

    size_t delay = (hilbert_taps - 1) / 2;
    phase_t phase = start_phase;
    for (size_t i = 0; i < n; ++i) {
        Q15 audio_delayed = (i >= delay) ? audio[i - delay] : Q15(0);
        Q15 hilbert_val = hilbert_out[i].re;
        Q15 re = kInvSqrt2 * audio_delayed;
        Q15 im = kInvSqrt2 * hilbert_val;
        Complex16 analytic = upper_sideband ? Complex16(re, im) : Complex16(re, -im);
        out[i] = analytic * sin_cos_q15(phase);
        phase = static_cast<phase_t>(phase + carrier_step);
    }
}

void ssb_demodulate(const Complex16* x, size_t n, phase_t start_phase, phase_t carrier_step, Q15* audio_out) {
    phase_t phase = start_phase;
    for (size_t i = 0; i < n; ++i) {
        Complex16 y = x[i] * sin_cos_q15(phase).conj();
        audio_out[i] = y.re;
        phase = static_cast<phase_t>(phase + carrier_step);
    }
}

} // namespace dsp
