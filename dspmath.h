#ifndef DSPMATH
#define DSPMATH

#include <dsptypes.h>

inline q15_t q15_mul(q15_t a, q15_t b);

inline q15_t q15_div(q15_t a, q15_t b);

inline complex16_t c16_add (complex16_t a, complex16_t b);

inline complex16_t c16_sub (complex16_t a, complex16_t b);

inline complex16_t c16_conj (complex16_t a, complex16_t b);

inline complex16_t c16_mul (complex16_t a, complex16_t b);

inline complex16_t c16_magsq (complex16_t a, complex16_t b);

inline complex16_t c16_phase (complex16_t a, complex16_t b);

inline q15_t sin_q15_cordic (q15_t phase);

inline q15_t cos_q15_cordic (q15_t phase);

inline q15_t* sin_q15_lut (size_t sz);

inline q15_t* cos_q15_lut (size_t sz);

#undef
