#include "../hilbert.hpp"
#include "../windows.hpp"

#include <cassert>
#include <cmath>
#include <vector>

namespace dsp {

void design_hilbert_q15(Q15* taps, size_t num_taps) {
    assert(num_taps > 0 && "design_hilbert_q15: num_taps must be > 0");
    assert(num_taps % 2 == 1 && "design_hilbert_q15: num_taps should be odd (Type III, zero center tap)");

    std::vector<Complex16> buf(num_taps);
    int center = static_cast<int>(num_taps / 2);
    for (size_t n = 0; n < num_taps; ++n) {
        int m = static_cast<int>(n) - center;
        double val = 0.0;
        if (m % 2 != 0) {
            val = 2.0 / (M_PI * static_cast<double>(m));
        }
        buf[n] = Complex16(Q15::from_double(val), Q15(0));
    }

    blackman_harris_window(buf.data(), num_taps, num_taps);

    for (size_t n = 0; n < num_taps; ++n) {
        taps[n] = buf[n].re;
    }
}

} // namespace dsp
