#ifndef DSP_TYPES_H
#define DSP_TYPES_H

#include <stdint.h>
#include <stddef.h>

typedef int16_t q15_t;

typedef struct complex16_t {
    q15_t re;
    q15_t im;
} complex16_t;

typedef uint16_t phase_t;

typedef struct sin_cos_q15_t {
    q15_t sin;
    q15_t cos;
} sin_cos_q15_t;

#define Q15_MAX ((q15_t)(0x7FFF))
#define Q15_MIN ((q15_t)(0x8000))

inline q15_t double_to_q15(double x); 

inline q15_t float_to_q15(float x); 
    
inline double q15_to_double(q15_t x);

inline float q15_to_float(q15_t x);

#endif
