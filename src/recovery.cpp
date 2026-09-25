#include "../recovery.hpp"
#include "../dspmath.hpp"

#include <cassert>

namespace dsp {

namespace {

Q15 signed_by(Q15 sign_of, Q15 value) {
    return (sign_of.raw() >= 0) ? value : -value;
}

} // namespace

GardnerTimingRecovery::GardnerTimingRecovery(Q15 kp, Q15 ki) : kp_(kp), ki_(ki) {}

bool GardnerTimingRecovery::step(Complex16 sample, Complex16* symbol_out) {
    hist_[0] = hist_[1];
    hist_[1] = sample;
    if (hist_count_ < 2) {
        ++hist_count_;
        return false;
    }

    Complex16 interp = hist_[0] + (hist_[1] - hist_[0]) * mu_;

    if (strobe_phase_ == 0) {
        last_half_ = interp;
        strobe_phase_ = 1;
        return false;
    }

    Complex16 on_time_new = interp;
    if (have_last_on_time_) {
        // Classic Gardner timing-error detector.
        Q15 error = ((on_time_new - last_on_time_) * last_half_.conj()).re;
        integrator_ = integrator_ + ki_ * error;
        Q15 adjustment = kp_ * error + integrator_;
        mu_ = mu_ - adjustment;
    }
    last_on_time_ = on_time_new;
    have_last_on_time_ = true;
    strobe_phase_ = 0;
    *symbol_out = on_time_new;
    return true;
}

CostasLoop::CostasLoop(int order, Q15 kp, Q15 ki) : order_(order), kp_(kp), ki_(ki) {
    assert(order == 2 || order == 4);
}

Complex16 CostasLoop::step(Complex16 sample) {
    Complex16 nco = sin_cos_q15(phase_);
    Complex16 y = sample * nco.conj();

    Q15 error;
    if (order_ == 2) {
        error = signed_by(y.re, y.im);
    } else {
        error = signed_by(y.re, y.im) - signed_by(y.im, y.re);
    }

    integrator_ = integrator_ + ki_ * error;
    Q15 adjustment = kp_ * error + integrator_;
    phase_ = static_cast<phase_t>(phase_ + static_cast<uint16_t>(adjustment.raw()));

    return y;
}

} // namespace dsp
