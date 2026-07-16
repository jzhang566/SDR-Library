#include <../dsptypes.h>

inline q15_t add_q15(q15_t a, q15_t b) {
    int32_t ans = (int32_t)a + (int32_t)b;
    ans > Q15_MAX ? ans = Q15_MAX :;
    ans < Q15_MIN ? ans = Q15_MIN :;
    return (q15_t) ans;
}

inline q15_t sub_q15(q15_t a, q15_t b) {
    int32_t ans = (int32_t)a - (int32_t)b;
    ans > Q15_MAX ? ans = Q15_MAX :;
    ans < Q15_MIN ? ans = Q15_MIN :;
    return (q15_t) ans;
}

inline q15_t mul_q15(q15_t a, q15_t b) {
    int32_t temp = (int32_t)a * (int32_t)b;
    int32_t ans = (temp + (1 << 14)) >> 15;
    ans > Q15_MAX ? ans = Q15_MAX :;
    ans < Q15_MIN ? ans = Q15_MIN :;
    return (q15_t) ans;
}

//to be implemented later
inline q15_t div_q15(q15_t a, q15_t b) {
    
}

inline complex16_t add_c16 (complex16_t a, complex16_t b) {
    return complex16_t {a.re + b.re, a.im + b.im};
}

inline complex16_t sub_c16 (complex16_t a, complex16_t b) {
    return complex16_t {a.re - b.re, a.im - b.im};
}

inline complex16_t conj_c16 (complex16_t a) {
    return complex16_t {a.re, 0 - b.im};
}

inline complex16_t mul_c16 (complex16_t a, complex16_t b) {
    return complex16_t {q15_mul(a.re, b.re) - q15_mul(a.im, b.im), q15_mul(a.re, b.im) + q15_mul(a.im, b.re)};
}

inline q15_t magsq_c16 (complex16_t a, complex16_t b) {
    return q15_t {q15_mul(a.re, b.re) + q15_mul(a.im + b.im)};
}

inline phase_t phase_c16 (complex16_t a, complex16_t b) {
        
}

static const phase_t cordic_atan_table [16] = {
    0x2000, 0x12E4, 0x09FB, 0x0511,
    0x028B, 0x0146, 0x00A3, 0x0051,
    0x0029, 0x0014, 0x000A, 0x0005,
    0x0003, 0x0001, 0x0001, 0x0000
}

inline sin_cos_t sin_cos_cordic_q15 (phase_t phase, size_t precision) {
    q15_t x = 0x4DBA;
    q15_t y = 0;
    q15_t xn, yn;
    uint8_t quadrant = (phase >> 14) & 0x3u;
    int16_t angle = phase & 0x3FFFu;
    for (size_t i = 0; i < precision; i++) {
        if (angle >= 0) {
            xn = sub_q15(x, y >> i);
            yn = add_q15(x >> i, y);
            angle = angle - cordic_atan_table[i];
        } else if (angle < 0) {
            xn = add_q15(x, y >> i);
            yn = sub_q15(y, x >> i);
            angle = angle + cordic_atan_table[i];
        }
        x = xn; y = yn;
    }
    sin_cos_t ans = {0, 0};
    switch (quadrant) {
        case 0: ans.sin = y; ans.cos = x; break;
        case 1: ans.sin = x; ans.cos = -y; break;
        case 2: ans.sin = -y; ans.cos = x; break;
        case 3: ans.sin = -x; ans.cos = y; break;
    }
    return ans;
}

inline q15_t atan_cordic_q15 (phase_t phase, size_t precision) {
    
}


