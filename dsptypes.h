#ifndef DSP_TYPES_H
#define DSP_TYPES_H

#include <stdint.h>
#include <stddef.h>

typedef int16_t dsp_Q15;
typedef int32_t dsp_Q31;

static inline dsp_Q15 double_to_Q15(double x) 

static inline dsp_Q31 double_to_Q31(double x) 

static inline dsp_Q15 float_to_Q15(float x) 

static inline dsp_Q31 float_to_Q31(float x) 
    
static inline double Q15_to_double(Q15 x)

static inline double Q31_to_double(Q31 x)

static inline float Q15_to_double(Q15 x)

static inline float Q31_to_double(Q31 x)

#endif
