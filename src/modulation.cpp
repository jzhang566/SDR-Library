#include "../modulation.hpp"

#include <cassert>

namespace dsp {

namespace {

constexpr Q15 kUnity(Q15::MAX); // symmetric "unit" amplitude (~0.99997)

uint32_t binary_to_gray(uint32_t b) {
    return b ^ (b >> 1);
}

uint32_t gray_to_binary(uint32_t g) {
    g ^= g >> 16;
    g ^= g >> 8;
    g ^= g >> 4;
    g ^= g >> 2;
    g ^= g >> 1;
    return g;
}

// Rounds num/den to the nearest integer, ties away from zero. den > 0.
int32_t round_div(int32_t num, int32_t den) {
    if (num >= 0) return (num + den / 2) / den;
    return -((-num + den / 2) / den);
}

// PAM levels per axis: QPSK is "4-QAM" (2 levels/axis); BPSK is the one
// 1-D scheme and isn't handled by this table.
unsigned levels_per_axis(ModScheme scheme) {
    switch (scheme) {
        case ModScheme::Qpsk:   return 2;
        case ModScheme::Qam16:  return 4;
        case ModScheme::Qam64:  return 8;
        case ModScheme::Qam256: return 16;
        case ModScheme::Bpsk:   return 0;
    }
    return 0;
}

unsigned log2_exact(unsigned v) {
    unsigned r = 0;
    while ((1u << r) < v) ++r;
    return r;
}

// Spacing between adjacent PAM levels' raw Q15 values. Bounds the corner
// points' VECTOR magnitude (not just their per-axis amplitude) to kUnity:
// corners sit at (+-(L-1)*step, +-(L-1)*step), vector magnitude
// sqrt(2)*(L-1)*step. Per-axis-only bounding would let a corner point's
// vector magnitude reach ~1.414x Q15 range, which is fine for static
// transmission but saturates the moment anything multiplies the symbol by
// a unit-magnitude rotator (e.g. Costas derotation can rotate that full
// vector magnitude onto a single axis). Uses 0.7 as a deliberately
// conservative approximation of 1/sqrt(2) (~0.7071) to keep this a simple,
// easily-verified-safe integer computation, at the cost of ~1% of the
// available dynamic range. L is a level *count* here, not a Q15 fraction,
// so this is plain integer arithmetic, not Q15::operator*.
int32_t axis_step(unsigned L) {
    int32_t per_axis = kUnity.raw() / static_cast<int32_t>(L - 1);
    return (per_axis * 7) / 10;
}

// index in [0, L) -> signed PAM amplitude, as a raw Q15 value.
int32_t axis_level(unsigned L, uint32_t index) {
    int32_t step = axis_step(L);
    return (2 * static_cast<int32_t>(index) - static_cast<int32_t>(L - 1)) * step;
}

// raw axis sample -> nearest level index in [0, L).
uint32_t axis_slice(unsigned L, int32_t raw) {
    int32_t step = axis_step(L);
    int32_t numerator = raw + step * static_cast<int32_t>(L - 1);
    int32_t index = round_div(numerator, 2 * step);
    if (index < 0) index = 0;
    if (index > static_cast<int32_t>(L - 1)) index = static_cast<int32_t>(L - 1);
    return static_cast<uint32_t>(index);
}

// First half of symbol_bits' bits -> I axis, second half -> Q axis; purely
// an internal convention (map_symbol/demap_symbol only need to agree with
// each other, not an external standard).
Complex16 map_square_qam(unsigned L, uint32_t symbol_bits) {
    unsigned bits_per_axis = log2_exact(L);
    uint32_t mask = (1u << bits_per_axis) - 1;
    uint32_t i_gray = (symbol_bits >> bits_per_axis) & mask;
    uint32_t q_gray = symbol_bits & mask;
    uint32_t i_index = gray_to_binary(i_gray);
    uint32_t q_index = gray_to_binary(q_gray);
    Q15 re(static_cast<int16_t>(axis_level(L, i_index)));
    Q15 im(static_cast<int16_t>(axis_level(L, q_index)));
    return Complex16(re, im);
}

uint32_t demap_square_qam(unsigned L, const Complex16& sample) {
    unsigned bits_per_axis = log2_exact(L);
    uint32_t i_index = axis_slice(L, sample.re.raw());
    uint32_t q_index = axis_slice(L, sample.im.raw());
    uint32_t i_gray = binary_to_gray(i_index);
    uint32_t q_gray = binary_to_gray(q_index);
    return (i_gray << bits_per_axis) | q_gray;
}

} // namespace

unsigned bits_per_symbol(ModScheme scheme) {
    switch (scheme) {
        case ModScheme::Bpsk:   return 1;
        case ModScheme::Qpsk:   return 2;
        case ModScheme::Qam16:  return 4;
        case ModScheme::Qam64:  return 6;
        case ModScheme::Qam256: return 8;
    }
    return 0;
}

Complex16 map_symbol(ModScheme scheme, uint32_t symbol_bits) {
    if (scheme == ModScheme::Bpsk) {
        return (symbol_bits & 0x1u) ? Complex16(-kUnity, Q15(0)) : Complex16(kUnity, Q15(0));
    }
    return map_square_qam(levels_per_axis(scheme), symbol_bits);
}

uint32_t demap_symbol(ModScheme scheme, Complex16 sample) {
    if (scheme == ModScheme::Bpsk) {
        return sample.re.raw() < 0 ? 1u : 0u;
    }
    return demap_square_qam(levels_per_axis(scheme), sample);
}

void modulate_bits(ModScheme scheme, const uint8_t* bits, size_t num_bits, Complex16* symbols_out) {
    unsigned bps = bits_per_symbol(scheme);
    assert(num_bits % bps == 0 && "modulate_bits: num_bits must be a multiple of bits_per_symbol(scheme)");
    size_t num_symbols = num_bits / bps;
    size_t bit_pos = 0; // MSB-first bit index into `bits`
    for (size_t s = 0; s < num_symbols; ++s) {
        uint32_t symbol_bits = 0;
        for (unsigned b = 0; b < bps; ++b) {
            size_t byte_idx = bit_pos / 8;
            unsigned bit_idx = 7 - static_cast<unsigned>(bit_pos % 8);
            uint32_t bit = (bits[byte_idx] >> bit_idx) & 0x1u;
            symbol_bits = (symbol_bits << 1) | bit;
            ++bit_pos;
        }
        symbols_out[s] = map_symbol(scheme, symbol_bits);
    }
}

size_t demodulate_bits(ModScheme scheme, const Complex16* symbols, size_t num_symbols, uint8_t* bits_out) {
    unsigned bps = bits_per_symbol(scheme);
    size_t total_bits = num_symbols * bps;
    size_t total_bytes = (total_bits + 7) / 8;
    for (size_t i = 0; i < total_bytes; ++i) bits_out[i] = 0;

    size_t bit_pos = 0;
    for (size_t s = 0; s < num_symbols; ++s) {
        uint32_t symbol_bits = demap_symbol(scheme, symbols[s]);
        for (unsigned b = 0; b < bps; ++b) {
            uint32_t bit = (symbol_bits >> (bps - 1 - b)) & 0x1u;
            size_t byte_idx = bit_pos / 8;
            unsigned bit_idx = 7 - static_cast<unsigned>(bit_pos % 8);
            bits_out[byte_idx] |= static_cast<uint8_t>(bit << bit_idx);
            ++bit_pos;
        }
    }
    return total_bits;
}

} // namespace dsp
