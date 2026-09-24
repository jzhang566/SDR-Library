#include "../fft.hpp"

#include <cassert>

namespace dsp {

void fft(uint8_t /*config*/, const Complex16* /*data*/, size_t /*N*/, Complex16* /*out*/, bool /*scale*/) {
    assert(false && "fft: not yet implemented");
}

} // namespace dsp
