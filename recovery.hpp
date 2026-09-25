#ifndef RECOVERY_HPP
#define RECOVERY_HPP

#include "dsptypes.hpp"

namespace dsp {

/* Gardner symbol-timing recovery, for a fixed 2 samples/symbol input.
 * Tracks a sub-sample interpolation offset (mu) via linear interpolation
 * between the two most recent input samples and the classic Gardner
 * timing-error detector, corrected each symbol by a PI loop filter.
 *
 * Scope limitation: mu is clamped to [0, ~1) rather than wrapping with a
 * true modulo-1 sample slip/repeat, so this converges for a FIXED
 * sub-sample timing offset (up to roughly +-0.5 samples) but does not
 * track a genuinely drifting sample-clock rate (frequency offset). A
 * polyphase interpolator + real slip logic would remove this limit;
 * left for a later pass.
 *
 * kp/ki are exposed directly rather than derived from a bandwidth/
 * damping-factor abstraction (the usual textbook loop-filter
 * parameterization): deriving Gardner-detector gain normalization
 * analytically risks a subtly wrong scale constant without empirical
 * calibration for this specific fixed-point detector, so tune these
 * directly instead. Smaller = slower but more stable lock. */
class GardnerTimingRecovery {
public:
    explicit GardnerTimingRecovery(Q15 kp = Q15::from_double(0.05), Q15 ki = Q15::from_double(0.0025));

    /* Feed one 2-sps input sample. Returns true and writes the recovered
     * (interpolated) on-time symbol to *symbol_out when a symbol strobe
     * fires this call (every other call, once locked). */
    bool step(Complex16 sample, Complex16* symbol_out);

    /* Current fractional interpolation offset, for monitoring/tuning
     * (e.g. plotting convergence). */
    Q15 mu() const { return mu_; }

private:
    Complex16 hist_[2];
    int hist_count_ = 0;
    Q15 mu_{0};
    int strobe_phase_ = 0;
    Complex16 last_half_;
    Complex16 last_on_time_;
    bool have_last_on_time_ = false;
    Q15 kp_, ki_;
    Q15 integrator_{0};
};

/* Costas loop carrier/phase recovery. order = 2 for BPSK (drives the
 * residual imaginary component to zero, decision-directed via the sign
 * of the real part), order = 4 for QPSK/QAM (the classic sign-sign
 * decision-directed detector). NCO uses sin_cos_q15() (honors the
 * CORDIC/LUT mode switch).
 *
 * Like GardnerTimingRecovery, kp/ki are exposed directly rather than via
 * a bandwidth/damping abstraction, for the same reason. The loop
 * filter's output is a Q15 value whose raw bits are reinterpreted
 * directly as a phase_t increment, so a saturated (+-1.0) correction is
 * a half rotation -- naturally well-scaled without a separate unit
 * conversion, and consistent with the wraparound phase-accumulator
 * pattern used elsewhere in this library (e.g. the FFT's twiddle
 * phases). */
class CostasLoop {
public:
    explicit CostasLoop(int order, Q15 kp = Q15::from_double(0.03), Q15 ki = Q15::from_double(0.001));

    /* Derotates one sample by the current phase estimate and updates the
     * loop from the phase-detector error. Returns the derotated sample. */
    Complex16 step(Complex16 sample);

    phase_t phase() const { return phase_; }

private:
    int order_;
    phase_t phase_ = 0;
    Q15 kp_, ki_;
    Q15 integrator_{0};
};

} // namespace dsp

#endif // RECOVERY_HPP
