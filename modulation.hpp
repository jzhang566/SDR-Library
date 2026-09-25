#ifndef MODULATION_HPP
#define MODULATION_HPP

#include "dsptypes.hpp"

namespace dsp {

enum class ModScheme {
    Bpsk,
    Qpsk,
    Qam16,
    Qam64,
    Qam256,
};

/* Number of bits carried per symbol for a scheme. */
unsigned bits_per_symbol(ModScheme scheme);

/* Maps the low bits_per_symbol(scheme) bits of symbol_bits to a Gray-coded
 * constellation point. Constellations are peak-normalized (the outermost
 * points sit at the largest representable Q15 magnitude per axis) rather
 * than unit-average-energy, so higher-order QAM's corner points never
 * saturate -- unit-average-energy 256-QAM would need corner points at
 * ~1.15x Q15 range. QPSK/16/64/256-QAM are all instances of one generic
 * square-constellation mapper (QPSK = "4-QAM", 2 levels/axis); BPSK is
 * the one 1-D special case. */
Complex16 map_symbol(ModScheme scheme, uint32_t symbol_bits);

/* Hard-decision slicer: nearest constellation point, returned as the
 * bits_per_symbol(scheme)-bit Gray code that produced it. O(1) per axis
 * (threshold arithmetic), not a search over the constellation. */
uint32_t demap_symbol(ModScheme scheme, Complex16 sample);

/* Maps a packed bit buffer (MSB-first within each byte) to symbols.
 * num_bits must be a multiple of bits_per_symbol(scheme). symbols_out
 * must hold num_bits / bits_per_symbol(scheme) elements. */
void modulate_bits(ModScheme scheme, const uint8_t* bits, size_t num_bits, Complex16* symbols_out);

/* Inverse of modulate_bits(): slices num_symbols symbols and packs their
 * Gray-coded bits (MSB-first) into bits_out, which must hold at least
 * ceil(num_symbols * bits_per_symbol(scheme) / 8) bytes. Returns the
 * number of bits written. */
size_t demodulate_bits(ModScheme scheme, const Complex16* symbols, size_t num_symbols, uint8_t* bits_out);

} // namespace dsp

#endif // MODULATION_HPP
