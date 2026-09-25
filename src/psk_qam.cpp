#include "../psk_qam.hpp"
#include "../filters.hpp"
#include "../pulse_shaping.hpp"

#include <cassert>

namespace dsp {

namespace {

std::vector<Complex16> make_rrc_taps(double rolloff, size_t sps, size_t rrc_taps) {
    std::vector<Q15> q15_taps(rrc_taps);
    design_rrc_q15(q15_taps.data(), rrc_taps, rolloff, sps);
    std::vector<Complex16> taps(rrc_taps);
    for (size_t i = 0; i < rrc_taps; ++i) taps[i] = Complex16(q15_taps[i], Q15(0));
    return taps;
}

} // namespace

PskQamModulator::PskQamModulator(ModScheme scheme, size_t sps, double rolloff, size_t rrc_taps)
    : scheme_(scheme), sps_(sps), rrc_taps_(make_rrc_taps(rolloff, sps, rrc_taps)) {}

size_t PskQamModulator::output_length(size_t num_bits) const {
    unsigned bps = bits_per_symbol(scheme_);
    assert(num_bits % bps == 0 && "PskQamModulator::output_length: num_bits must be a multiple of bits_per_symbol(scheme)");
    return (num_bits / bps) * sps_;
}

void PskQamModulator::modulate(const uint8_t* bits, size_t num_bits, Complex16* out) {
    unsigned bps = bits_per_symbol(scheme_);
    assert(num_bits % bps == 0 && "PskQamModulator::modulate: num_bits must be a multiple of bits_per_symbol(scheme)");
    size_t num_symbols = num_bits / bps;

    std::vector<Complex16> symbols(num_symbols);
    modulate_bits(scheme_, bits, num_bits, symbols.data());

    std::vector<Complex16> upsampled(num_symbols * sps_);
    interpolate(symbols.data(), upsampled.data(), sps_, num_symbols);

    fir_filter(upsampled.data(), rrc_taps_.data(), out, upsampled.size(), rrc_taps_.size());
}

PskQamDemodulator::PskQamDemodulator(ModScheme scheme, size_t sps, double rolloff, size_t rrc_taps,
                                      Q15 gardner_kp, Q15 gardner_ki, Q15 costas_kp, Q15 costas_ki)
    : scheme_(scheme),
      sps_(sps),
      rrc_taps_(make_rrc_taps(rolloff, sps, rrc_taps)),
      gardner_(gardner_kp, gardner_ki),
      costas_(scheme == ModScheme::Bpsk ? 2 : 4, costas_kp, costas_ki) {}

size_t PskQamDemodulator::max_output_bits(size_t n) const {
    return (n / sps_ + 1) * bits_per_symbol(scheme_);
}

size_t PskQamDemodulator::demodulate(const Complex16* x, size_t n, uint8_t* bits_out) {
    std::vector<Complex16> filtered(n);
    fir_filter(x, rrc_taps_.data(), filtered.data(), n, rrc_taps_.size());

    std::vector<Complex16> symbols;
    symbols.reserve(n / sps_ + 1);
    for (size_t i = 0; i < n; ++i) {
        Complex16 sym;
        if (gardner_.step(filtered[i], &sym)) {
            symbols.push_back(costas_.step(sym));
        }
    }

    return demodulate_bits(scheme_, symbols.data(), symbols.size(), bits_out);
}

} // namespace dsp
