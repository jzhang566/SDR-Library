#include "../windows.hpp"
#include "../dspmath.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace dsp {

namespace {

// Maps n in [0, sz) to a phase_t angle covering [0, 2*pi) as n/(sz-1)
// goes from 0 to 1, i.e. the phase_t units for 2*pi*n/(sz-1).
phase_t window_phase(size_t n, size_t sz) {
    return static_cast<phase_t>((static_cast<uint64_t>(n) * 0x10000u) / (sz - 1));
}

void zero_pad(Complex16* x, size_t sz, size_t N) {
    for (size_t n = sz; n < N; ++n) x[n] = Complex16(Q15(0), Q15(0));
}

// I0(x) = sum_{k=0}^inf ((x/2)^(2k)) / (k!)^2, evaluated in double since
// I0(beta) for typical Kaiser beta values (0-20ish) is far outside Q15
// range; only the final normalized ratio gets quantized to Q15.
double bessel_i0(double x) {
    double sum = 1.0;
    double term = 1.0;
    double half_x_sq = (x * 0.5) * (x * 0.5);
    for (int k = 1; k < 32; ++k) {
        term *= half_x_sq / (static_cast<double>(k) * static_cast<double>(k));
        sum += term;
        if (term < 1e-15 * sum) break;
    }
    return sum;
}

} // namespace

void hamming_window(Complex16* x, size_t sz, size_t N) {
    assert(sz > 0 && sz <= N && "hamming_window: need 0 < sz <= N");
    constexpr Q15 a0(17695); // round(0.54   * 32768)
    constexpr Q15 a1(15073); // round(0.46   * 32768)
    if (sz > 1) {
        for (size_t n = 0; n < sz; ++n) {
            Q15 c = sin_cos_q15(window_phase(n, sz)).re;
            Q15 w = a0 - a1 * c;
            x[n] = x[n] * w;
        }
    }
    zero_pad(x, sz, N);
}

void blackman_harris_window(Complex16* x, size_t sz, size_t N) {
    assert(sz > 0 && sz <= N && "blackman_harris_window: need 0 < sz <= N");
    constexpr Q15 a0(11756); // round(0.35875 * 32768)
    constexpr Q15 a1(16000); // round(0.48829 * 32768)
    constexpr Q15 a2(4629);  // round(0.14128 * 32768)
    constexpr Q15 a3(383);   // round(0.01168 * 32768)
    if (sz > 1) {
        for (size_t n = 0; n < sz; ++n) {
            phase_t p1 = window_phase(n, sz);
            phase_t p2 = static_cast<phase_t>(p1 * 2);
            phase_t p3 = static_cast<phase_t>(p1 * 3);
            Q15 c1 = sin_cos_q15(p1).re;
            Q15 c2 = sin_cos_q15(p2).re;
            Q15 c3 = sin_cos_q15(p3).re;
            Q15 w = a0 - a1 * c1 + a2 * c2 - a3 * c3;
            x[n] = x[n] * w;
        }
    }
    zero_pad(x, sz, N);
}

void kaiser_window(Complex16* x, size_t sz, size_t N, double beta) {
    assert(sz > 0 && sz <= N && "kaiser_window: need 0 < sz <= N");
    if (sz > 1) {
        double i0_beta = bessel_i0(beta);
        double denom = static_cast<double>(sz - 1);
        for (size_t n = 0; n < sz; ++n) {
            double r = (2.0 * static_cast<double>(n) / denom) - 1.0;
            double arg = beta * std::sqrt(std::max(0.0, 1.0 - r * r));
            Q15 w = Q15::from_double(bessel_i0(arg) / i0_beta);
            x[n] = x[n] * w;
        }
    }
    zero_pad(x, sz, N);
}

} // namespace dsp
