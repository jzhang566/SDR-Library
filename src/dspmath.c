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
    return complex16_t {a.re, 0 - a.im};
}

inline complex16_t mul_c16 (complex16_t a, complex16_t b) {
    return complex16_t {add_q15(
        sub_q15(mul_q15(a.re, b.re), mul_q15(a.im, b.im)), 
        add_q15(mul_q15(a.re, b.im), mul_q15(a.im, b.re))};
}

inline q15_t magsq_c16 (complex16_t a, complex16_t b) {
    return q15_t {add_q15(mul_q15(a.re, b.re), mul_q15(a.im + b.im))};
}

inline phase_t phase_c16 (complex16_t a) {
    return atan_cordic_q15(a);
}

static const phase_t cordic_atan_table [16] = {
    0x2000, 0x12E4, 0x09FB, 0x0511,
    0x028B, 0x0146, 0x00A3, 0x0051,
    0x0029, 0x0014, 0x000A, 0x0005,
    0x0003, 0x0001, 0x0001, 0x0000
}

static int cordic_precision = 16;

inline sin_cos_t sin_cos_cordic_q15 (phase_t phase) {
    q15_t x = 0x4DBA;
    q15_t y = 0;
    q15_t xn, yn;
    uint8_t quadrant = (phase >> 14) & 0x3u;
    int16_t angle = phase & 0x3FFFu;
    for (size_t i = 0; i < cordic_precision; i++) {
        if (angle >= 0) {
            xn = sub_q15(x, y >> i);
            yn = add_q15(y, x >> i);
            angle = angle - cordic_atan_table[i];
        } else {
            xn = add_q15(x, y >> i);
            yn = sub_q15(y, x >> i);
            angle = angle + cordic_atan_table[i];
        }
        x = xn; y = yn;
    }
    sin_cos_t ans = {0, 0};
    switch (quadrant) {
        case 0: return sin_cos_t {y, x};
        case 1: return sin_cos_t {x, -y};
        case 2: return sin_cos_t {-y, x};
        case 3: return sin_cos_t {-x, y};
    }
}

inline phase_t atan_cordic_q15 (complex16_t v) {
    q15_t x = v.re;
    q15_t y = v.im;
    phase_t angle = 0;
    if (x < 0) {
        x = sub_q15(0, x);
        y = sub_q15(0, y);
        angle = 0x8000;
    }
    q15_t xn, yn;
    for (size_t i = 0; i < cordic_precision; i++) {
        if (y >= 0) {
            xn = add_q15(x, y >> i);
            yn = sub_q15(y, x >> i);
            angle = angle + cordic_atan_table[i];
        } else {
            xn = sub_q15(x, y >> i);
            yn = add_q15(y, x >> i);
            angle = angle - cordic_atan_table[i];
        }
        x = xn;
        y = yn;
    }
    return angle; 
}

typedef struct sinusoid_lut_q15_t {
    q15_t *lut;
    uint16_t sz;
    uint16_t frac_bits;
    uint8_t bits;
    
    
} sinusoid_lut_q15_t;
static sinusoid_lut_q15_t cos_lut_q15_g;
// ptr is a user-allocated buffer of size (sz * sizeof(q15_t))
// Size must be a power of 2, greater than 4. This allows for quick generation and lookup without division.
// Only one sin-cos LUT can be created at a time. 
inline bool generate_sinusoid_lut_q15 (q15_t *ptr, size_t sz) {
    if (sz < 8 || (sz & (sz - 1)) != 0) return 1;
    for (size_t i = 0; i < sz; i++) {
        if ((sz >> i) & 0x1) {
            cos_lut_q15_g.bits = i;
            break;
        }
    }
    cos_lut_q15_g.sz = sz;
    cos_lut_q15_g.lut = ptr;
    cos_lut_q15_g.frac_bits = 14 - cos_lut_q15_g.bits;
    size_t ind = 0;
    for (size_t i = 0; i < 0x4000; i += 1u << cos_lut_q15_g.frac_bits) {
        ptr[ind] = sin_cos_cordic_q15(i).cos;
        ind = ind + 1;
    }
    return 0;
}

inline q15_t cos_lut_q15 (phase_t angle) {
    uint8_t quadrant = (angle >> 14) & 0x3u;
    uint16_t r = angle & 0x3FFFu;
    bool neg = 0;
    switch (quadrant) {
        case 0: break;
        case 1: neg = 1; r = (0x4000u - r); break;
        case 2: neg = 1; break;
        case 3: r = (0x4000u - r); break;
    }
    uint16_t ind = r >> cos_lut_q15_g.frac_bits;
    uint16_t frac = r & ((1u << (cos_lut_q15_g.frac_bits)) - 1);
    q15_t frac_q15 = (uint32_t)(frac << (15 - cos_lut_q15_g.frac_bits));
    q15_t y0 = cos_lut_q15_g.lut[ind];
    q15_t y1 = (ind == cos_lut_q15_g.sz - 1) ? 0 : cos_lut_q15_g.lut[ind + 1];
    q15_t diff = sub_q15(y1, y0);
    q15_t interp = mul_q15(diff, frac_q15);
    return neg ? -1 * add_q15(y0, interp) : add_q15(y0, interp);
}

inline q15_t sin_lut_q15 (phase_t angle) {
    return cos_lut_q15((0x4000u - phase));
}

inline q15_t cos_lut_q15_nointerp (phase_t angle) {
    uint8_t quadrant = (angle >> 14) & 0x3u;
    uint16_t r = angle & 0x3FFFu;
    bool neg = 0;
    switch (quadrant) {
        case 0: break;
        case 1: neg = 1; r = (0x4000u - r); break;
        case 2: neg = 1; break;
        case 3: r = (0x4000u - r); break;
    }

    }
    uint16_t ind = (r + (1 << (cos_lut_q15_g.frac_bits - 1))) >> cos_lut_q15_g.frac_bits;
    return neg ? -1 * cos_lut_q15_g.lut[ind] : cos_lut_q15_g.lut[ind];
}

inline q15_t sin_lut_q15_nointerp (phase_t angle) {
    return cos_lut_q15_nointerp (angle);
}
