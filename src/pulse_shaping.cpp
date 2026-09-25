#include "../pulse_shaping.hpp"

#include <cassert>
#include <cmath>
#include <vector>

namespace dsp {

void design_rrc_q15(Q15* taps, size_t num_taps, double rolloff, size_t sps) {
    assert(num_taps > 0 && "design_rrc_q15: num_taps must be > 0");
    assert(num_taps % 2 == 1 && "design_rrc_q15: num_taps should be odd for a symmetric filter");
    assert(sps > 0 && "design_rrc_q15: sps must be > 0");
    assert(rolloff >= 0.0 && rolloff <= 1.0 && "design_rrc_q15: rolloff must be in [0, 1]");

    const double beta = rolloff;
    const double center = (static_cast<double>(num_taps) - 1.0) / 2.0;

    std::vector<double> h(num_taps);
    for (size_t n = 0; n < num_taps; ++n) {
        double x = (static_cast<double>(n) - center) / static_cast<double>(sps); // t/T

        double val;
        if (beta > 0.0 && std::fabs(std::fabs(x) - 1.0 / (4.0 * beta)) < 1e-9) {
            // Analytic limit at t = +-T/(4*beta), the RRC formula's other
            // removable singularity (0/0 there for a raw evaluation).
            val = (beta / std::sqrt(2.0)) *
                  ((1.0 + 2.0 / M_PI) * std::sin(M_PI / (4.0 * beta)) +
                   (1.0 - 2.0 / M_PI) * std::cos(M_PI / (4.0 * beta)));
        } else if (std::fabs(x) < 1e-9) {
            // Analytic limit at t = 0.
            val = 1.0 - beta + 4.0 * beta / M_PI;
        } else {
            double num = std::sin(M_PI * x * (1.0 - beta)) + 4.0 * beta * x * std::cos(M_PI * x * (1.0 + beta));
            double den = M_PI * x * (1.0 - (4.0 * beta * x) * (4.0 * beta * x));
            val = num / den;
        }
        h[n] = val;
    }

    double peak = 0.0;
    for (double v : h) peak = std::max(peak, std::fabs(v));

    double target = static_cast<double>(Q15::MAX) / 32768.0;
    double scale = (peak > 0.0) ? (target / peak) : 1.0;
    for (size_t n = 0; n < num_taps; ++n) {
        taps[n] = Q15::from_double(h[n] * scale);
    }
}

} // namespace dsp
