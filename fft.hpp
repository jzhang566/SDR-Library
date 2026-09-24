#ifndef FFT_HPP
#define FFT_HPP

#include "dsptypes.hpp"

namespace dsp {

/* Not yet implemented. config/scale semantics (radix, in-place vs.
 * out-of-place, forward/inverse) are still undecided; calling this aborts. */
void fft(uint8_t config, const Complex16* data, size_t N, Complex16* out, bool scale);

} // namespace dsp

#endif // FFT_HPP
