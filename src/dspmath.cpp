#include "../dspmath.hpp"

#include <cassert>

namespace dsp {

namespace {

constexpr int kCordicPrecision = 16;

constexpr phase_t kCordicAtanTable[16] = {
    0x2000, 0x12E4, 0x09FB, 0x0511,
    0x028B, 0x0146, 0x00A3, 0x0051,
    0x0029, 0x0014, 0x000A, 0x0005,
    0x0003, 0x0001, 0x0001, 0x0000
};

struct SinusoidLutState {
    Q15* lut = nullptr;
    uint16_t sz = 0;
    uint16_t frac_bits = 0;
    uint8_t bits = 0;
};

SinusoidLutState g_lut;

SinCosMode g_sin_cos_mode = SinCosMode::Cordic;

} // namespace

phase_t phase_c16(const Complex16& a) {
    return atan_cordic_q15(a);
}

Complex16 sin_cos_cordic_q15(phase_t phase) {
    Q15 x(0x4DBA);
    Q15 y(0);
    Q15 xn, yn;
    uint8_t quadrant = (phase >> 14) & 0x3u;
    int16_t angle = static_cast<int16_t>(phase & 0x3FFFu);
    for (int i = 0; i < kCordicPrecision; ++i) {
        if (angle >= 0) {
            xn = x - (y >> i);
            yn = y + (x >> i);
            angle = static_cast<int16_t>(angle - kCordicAtanTable[i]);
        } else {
            xn = x + (y >> i);
            yn = y - (x >> i);
            angle = static_cast<int16_t>(angle + kCordicAtanTable[i]);
        }
        x = xn;
        y = yn;
    }
    switch (quadrant) {
        case 0: return Complex16{x, y};
        case 1: return Complex16{-y, x};
        case 2: return Complex16{-x, -y};
        default: return Complex16{y, -x};
    }
}

phase_t atan_cordic_q15(const Complex16& v) {
    Q15 x = v.re;
    Q15 y = v.im;
    phase_t angle = 0;
    if (x < Q15(0)) {
        x = -x;
        y = -y;
        angle = 0x8000;
    }
    Q15 xn, yn;
    for (int i = 0; i < kCordicPrecision; ++i) {
        if (y >= Q15(0)) {
            xn = x + (y >> i);
            yn = y - (x >> i);
            angle = static_cast<phase_t>(angle + kCordicAtanTable[i]);
        } else {
            xn = x - (y >> i);
            yn = y + (x >> i);
            angle = static_cast<phase_t>(angle - kCordicAtanTable[i]);
        }
        x = xn;
        y = yn;
    }
    return angle;
}

// sz must be a power of 2, greater than 4. This allows for quick generation
// and lookup without division.
bool generate_sinusoid_lut_q15(Q15* ptr, size_t sz) {
    if (sz < 8 || (sz & (sz - 1)) != 0) return false;

    uint8_t bits = 0;
    for (size_t i = 0; i < sz; ++i) {
        if ((sz >> i) & 0x1u) {
            bits = static_cast<uint8_t>(i);
            break;
        }
    }

    g_lut.lut = ptr;
    g_lut.sz = static_cast<uint16_t>(sz);
    g_lut.bits = bits;
    g_lut.frac_bits = static_cast<uint16_t>(14 - bits);

    uint16_t ind = 0;
    for (uint16_t i = 0; i < 0x4000; i = static_cast<uint16_t>(i + (1u << g_lut.frac_bits))) {
        ptr[ind] = sin_cos_cordic_q15(i).re;
        ++ind;
    }
    return true;
}

Q15 cos_lut_q15(phase_t angle) {
    uint8_t quadrant = (angle >> 14) & 0x3u;
    uint16_t r = angle & 0x3FFFu;
    bool neg = false;
    switch (quadrant) {
        case 0: break;
        case 1: neg = true; r = static_cast<uint16_t>(0x4000u - r); break;
        case 2: neg = true; break;
        case 3: r = static_cast<uint16_t>(0x4000u - r); break;
    }
    uint16_t ind = static_cast<uint16_t>(r >> g_lut.frac_bits);
    uint16_t frac = static_cast<uint16_t>(r & ((1u << g_lut.frac_bits) - 1));
    Q15 frac_q15(static_cast<int16_t>(frac << (15 - g_lut.frac_bits)));
    Q15 y0 = g_lut.lut[ind];
    Q15 y1 = (ind == g_lut.sz - 1) ? Q15(0) : g_lut.lut[ind + 1];
    Q15 diff = y1 - y0;
    Q15 interp = diff * frac_q15;
    Q15 result = y0 + interp;
    return neg ? -result : result;
}

Q15 sin_lut_q15(phase_t angle) {
    return cos_lut_q15(static_cast<phase_t>(0x4000u - angle));
}

Q15 cos_lut_q15_nointerp(phase_t angle) {
    uint8_t quadrant = (angle >> 14) & 0x3u;
    uint16_t r = angle & 0x3FFFu;
    bool neg = false;
    switch (quadrant) {
        case 0: break;
        case 1: neg = true; r = static_cast<uint16_t>(0x4000u - r); break;
        case 2: neg = true; break;
        case 3: r = static_cast<uint16_t>(0x4000u - r); break;
    }
    uint16_t ind = static_cast<uint16_t>((r + (1u << (g_lut.frac_bits - 1))) >> g_lut.frac_bits);
    Q15 value = g_lut.lut[ind];
    return neg ? -value : value;
}

Q15 sin_lut_q15_nointerp(phase_t angle) {
    return cos_lut_q15_nointerp(angle);
}

void set_sin_cos_mode(SinCosMode mode) {
    g_sin_cos_mode = mode;
}

SinCosMode sin_cos_mode() {
    return g_sin_cos_mode;
}

Complex16 sin_cos_q15(phase_t phase) {
    switch (g_sin_cos_mode) {
        case SinCosMode::Lut:
            assert(g_lut.lut != nullptr &&
                   "sin_cos_q15: Lut mode selected but no LUT generated (call generate_sinusoid_lut_q15 first)");
            return Complex16(cos_lut_q15(phase), sin_lut_q15(phase));
        case SinCosMode::Cordic:
        default:
            return sin_cos_cordic_q15(phase);
    }
}

namespace {

// Integer sqrt via the classic digit-by-digit (binary search) method.
uint32_t isqrt32(uint32_t n) {
    uint32_t res = 0;
    uint32_t bit = 1u << 30; // highest power of 4 <= UINT32_MAX
    while (bit > n) bit >>= 2;
    while (bit != 0) {
        uint32_t trial = res + bit;
        if (n >= trial) {
            n -= trial;
            res = (res >> 1) + bit;
        } else {
            res >>= 1;
        }
        bit >>= 2;
    }
    return res;
}

} // namespace

Q15 sqrt_q15(Q15 x) {
    assert(x.raw() >= 0 && "sqrt_q15: x must be non-negative");
    // sqrt(v) * 32768 == isqrt(raw * 32768) for v = raw/32768, raw >= 0.
    uint32_t v = static_cast<uint32_t>(x.raw()) << 15;
    uint32_t r = isqrt32(v);
    if (r > static_cast<uint32_t>(Q15::MAX)) r = Q15::MAX;
    return Q15(static_cast<int16_t>(r));
}

void pointwise_add(const Complex16* a, const Complex16* b, Complex16* out, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = a[i] + b[i];
}

void pointwise_sub(const Complex16* a, const Complex16* b, Complex16* out, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = a[i] - b[i];
}

void pointwise_mul(const Complex16* a, const Complex16* b, Complex16* out, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = a[i] * b[i];
}

void pointwise_mul(const Complex16* a, const Q15* b, Complex16* out, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = a[i] * b[i];
}

void scale(const Complex16* a, Q15 factor, Complex16* out, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = a[i] * factor;
}

void magsq(const Complex16* a, Q15* out, size_t n) {
    for (size_t i = 0; i < n; ++i) out[i] = a[i].magsq();
}

namespace {

// Sums a[i] * conj(b[i]) in 64-bit fixed point (each raw product is Q1.30;
// the accumulator holds their unsaturated sum) and rounds/saturates once to
// Q15 at the end, instead of saturating after every term.
Complex16 dot_wide(const Complex16* a, const Complex16* b, size_t n) {
    int64_t re_acc = 0;
    int64_t im_acc = 0;
    for (size_t i = 0; i < n; ++i) {
        int64_t ar = a[i].re.raw();
        int64_t ai = a[i].im.raw();
        int64_t br = b[i].re.raw();
        int64_t bi = b[i].im.raw();
        re_acc += ar * br + ai * bi;
        im_acc += ai * br - ar * bi;
    }
    auto to_q15 = [](int64_t acc) -> Q15 {
        int64_t rounded = (acc + (1LL << 14)) >> 15;
        if (rounded > Q15::MAX) rounded = Q15::MAX;
        if (rounded < Q15::MIN) rounded = Q15::MIN;
        return Q15(static_cast<int16_t>(rounded));
    };
    return Complex16(to_q15(re_acc), to_q15(im_acc));
}

} // namespace

Complex16 dot(const Complex16* a, const Complex16* b, size_t n) {
    return dot_wide(a, b, n);
}

Q15 energy(const Complex16* a, size_t n) {
    // dot(a, a)'s imaginary part is exactly 0 (ai*ar - ar*ai), so this is
    // just sum_i |a[i]|^2.
    return dot_wide(a, a, n).re;
}

} // namespace dsp
