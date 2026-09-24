#include "../filters.hpp"

namespace dsp {

void fir_filter(const Complex16* x, const Complex16* h, Complex16* y, size_t sz, size_t num_taps) {
    for (size_t n = 0; n < sz; ++n) {
        Complex16 acc;
        size_t taps = (num_taps < n + 1) ? num_taps : n + 1;
        for (size_t k = 0; k < taps; ++k) {
            acc += h[k] * x[n - k];
        }
        y[n] = acc;
    }
}

void biquad_filter(const Complex16* x, Q15 b0, Q15 b1, Q15 b2, Q15 a1, Q15 a2, Complex16* y, size_t sz) {
    for (size_t t = 2; t < sz; ++t) {
        y[t] = b0 * x[t] + b1 * x[t - 1] + b2 * x[t - 2] - a1 * y[t - 1] - a2 * y[t - 2];
    }
}

} // namespace dsp
