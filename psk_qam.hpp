#ifndef PSK_QAM_HPP
#define PSK_QAM_HPP

#include "dsptypes.hpp"
#include "modulation.hpp"
#include "recovery.hpp"

#include <vector>

namespace dsp {

/* Ties together modulation.hpp's symbol mapper and pulse_shaping.hpp's
 * RRC filter into a usable TX chain: bits -> symbols -> zero-stuff
 * upsample -> RRC pulse shape. rrc_taps should be odd (see
 * design_rrc_q15()). */
class PskQamModulator {
public:
    PskQamModulator(ModScheme scheme, size_t sps, double rolloff, size_t rrc_taps);

    /* num_bits must be a multiple of bits_per_symbol(scheme). out must
     * hold output_length(num_bits) elements. */
    void modulate(const uint8_t* bits, size_t num_bits, Complex16* out);

    size_t output_length(size_t num_bits) const;

private:
    ModScheme scheme_;
    size_t sps_;
    std::vector<Complex16> rrc_taps_;
};

/* Ties together the RX matched RRC filter, GardnerTimingRecovery, and
 * CostasLoop into a usable RX chain: oversampled IQ -> matched filter ->
 * Gardner timing recovery (decimates to symbol rate) -> Costas carrier
 * recovery (runs at symbol rate, where a decision-directed detector is
 * sharpest) -> hard-decision slicing -> packed bits.
 *
 * Timing recovery runs before carrier recovery deliberately: Gardner's
 * detector isn't decision-directed and tolerates modest residual carrier
 * phase, whereas feeding Costas's actively-converging phase estimate
 * into Gardner first would corrupt the timing detector while the carrier
 * loop is still settling.
 *
 * Burst/block model (see the modulation-subsystem plan): each
 * demodulate() call is matched-filtered independently, so the first few
 * taps' worth of samples at the start of each call have reduced
 * accuracy (fir_filter's "samples before x[0] are absent" contract). The
 * Gardner/Costas loop state *does* persist across calls, so a long
 * capture fed in chunks stays locked. */
class PskQamDemodulator {
public:
    PskQamDemodulator(ModScheme scheme, size_t sps, double rolloff, size_t rrc_taps,
                       Q15 gardner_kp = Q15::from_double(0.05), Q15 gardner_ki = Q15::from_double(0.0025),
                       Q15 costas_kp = Q15::from_double(0.03), Q15 costas_ki = Q15::from_double(0.001));

    /* Feeds n oversampled input samples through the chain above. bits_out
     * must hold at least max_output_bits(n) / 8 (rounded up) bytes.
     * Returns the number of bits actually written.
     *
     * Total pipeline latency from a transmitted symbol to its recovered
     * bits is (rrc_taps - 1) / sps - 1 symbols: one symbol *less* than
     * the RRC matched filter's own group delay, since
     * GardnerTimingRecovery's strobe/interpolation bookkeeping ends up
     * sampling slightly ahead of that naive estimate. Confirmed
     * empirically via an end-to-end test with a known bit sequence --
     * verify this offset yourself if you change sps or the recovery
     * loops' internals, since it isn't purely a closed-form filter-delay
     * calculation. */
    size_t demodulate(const Complex16* x, size_t n, uint8_t* bits_out);

    size_t max_output_bits(size_t n) const;

private:
    ModScheme scheme_;
    size_t sps_;
    std::vector<Complex16> rrc_taps_;
    GardnerTimingRecovery gardner_;
    CostasLoop costas_;
};

} // namespace dsp

#endif // PSK_QAM_HPP
