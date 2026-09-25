#include "../fft.hpp"
#include "../dspmath.hpp"

#include <cassert>
#include <vector>

namespace dsp {

namespace {

// twiddles[k] = e^{-j*2*pi*k/N_top} (or +j... for inv) for k in
// [0, N_top/2), computed once by fft() for the top-level N. At recursion
// depth d (sub-problem size N_top/2^d), the twiddle needed for index k is
// twiddles[k * 2^d]: e^{-j*2*pi*k/(N_top/2^d)} == e^{-j*2*pi*(k*2^d)/N_top}.
// So every level reuses a strided slice of the same table instead of
// recomputing its own -- O(N) sin/cos calls total instead of O(N log N).
void fft_impl(const Complex16* x, size_t N, Complex16* X, char scale_mode,
              const Complex16* twiddles, size_t stride) {
    if (N == 1) {
        X[0] = x[0];
        return;
    }

    size_t N_half = N >> 1;
    std::vector<Complex16> even(N_half), odd(N_half);
    for (size_t i = 0; i < N_half; ++i) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    std::vector<Complex16> E(N_half), O(N_half);
    fft_impl(even.data(), N_half, E.data(), scale_mode, twiddles, stride * 2);
    fft_impl(odd.data(), N_half, O.data(), scale_mode, twiddles, stride * 2);

    for (size_t k = 0; k < N_half; ++k) {
        Complex16 t = O[k] * twiddles[k * stride];
        Complex16 sum = E[k] + t;
        Complex16 diff = E[k] - t;
        if (scale_mode == 's') {
            sum = sum >> 1;
            diff = diff >> 1;
        }
        X[k] = sum;
        X[k + N_half] = diff;
    }
}

} // namespace

void fft(Complex16* x, size_t N, Complex16* X, char scale_mode, bool inv) {
    assert(N > 0 && N <= 0x10000 && (N & (N - 1)) == 0 && "fft: N must be a power of two, <= 0x10000");
    assert((scale_mode == 's' || scale_mode == 'u') && "fft: unknown scale_mode");
    assert(scale_mode != 'u' && "fft: scale_mode 'u' (unscaled) not yet implemented");

    if (N == 1) {
        X[0] = x[0];
        return;
    }

    size_t N_half = N >> 1;
    std::vector<Complex16> twiddles(N_half);
    uint16_t step = static_cast<uint16_t>(0x10000u / N);
    for (size_t k = 0; k < N_half; ++k) {
        uint32_t idx = static_cast<uint32_t>(k) * step;
        phase_t phase = inv ? static_cast<phase_t>(idx) : static_cast<phase_t>(-idx);
        twiddles[k] = sin_cos_q15(phase);
    }

    fft_impl(x, N, X, scale_mode, twiddles.data(), 1);
}

} // namespace dsp
