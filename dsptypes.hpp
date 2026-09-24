#ifndef DSP_TYPES_HPP
#define DSP_TYPES_HPP

#include <cstdint>
#include <cstddef>

namespace dsp {

// Q1.15 fixed-point value with saturating arithmetic.
class Q15 {
public:
    static constexpr int16_t MAX = 0x7FFF;
    static constexpr int16_t MIN = -32768;

    constexpr Q15() noexcept : raw_(0) {}
    constexpr Q15(int16_t raw) noexcept : raw_(raw) {}

    static Q15 from_double(double x);
    static Q15 from_float(float x);

    constexpr operator int16_t() const noexcept { return raw_; }
    constexpr int16_t raw() const noexcept { return raw_; }
    double to_double() const;
    float to_float() const;

    Q15 operator+(Q15 rhs) const;
    Q15 operator-(Q15 rhs) const;
    Q15 operator*(Q15 rhs) const;
    Q15 operator/(Q15 rhs) const;
    Q15 operator-() const;
    Q15 operator>>(int shift) const noexcept { return Q15(static_cast<int16_t>(raw_ >> shift)); }
    Q15 operator<<(int shift) const noexcept { return Q15(static_cast<int16_t>(raw_ << shift)); }

    Q15& operator+=(Q15 rhs) { return *this = *this + rhs; }
    Q15& operator-=(Q15 rhs) { return *this = *this - rhs; }
    Q15& operator*=(Q15 rhs) { return *this = *this * rhs; }
    Q15& operator/=(Q15 rhs) { return *this = *this / rhs; }

    constexpr bool operator==(Q15 rhs) const noexcept { return raw_ == rhs.raw_; }
    constexpr bool operator!=(Q15 rhs) const noexcept { return raw_ != rhs.raw_; }
    constexpr bool operator<(Q15 rhs) const noexcept { return raw_ < rhs.raw_; }
    constexpr bool operator>(Q15 rhs) const noexcept { return raw_ > rhs.raw_; }
    constexpr bool operator<=(Q15 rhs) const noexcept { return raw_ <= rhs.raw_; }
    constexpr bool operator>=(Q15 rhs) const noexcept { return raw_ >= rhs.raw_; }

private:
    int16_t raw_;
};

using phase_t = uint16_t;

// Complex value with Q15 real/imaginary components, saturating arithmetic.
// Also used as a unit-magnitude phasor (re = cos, im = sin) per Euler's
// formula, e.g. as returned by sin_cos_cordic_q15 — multiplying a sample by
// one rotates its phase.
class Complex16 {
public:
    Q15 re;
    Q15 im;

    constexpr Complex16() noexcept : re(), im() {}
    constexpr Complex16(Q15 re_, Q15 im_) noexcept : re(re_), im(im_) {}

    Complex16 operator+(const Complex16& rhs) const { return {re + rhs.re, im + rhs.im}; }
    Complex16 operator-(const Complex16& rhs) const { return {re - rhs.re, im - rhs.im}; }
    Complex16 operator*(const Complex16& rhs) const {
        return {re * rhs.re - im * rhs.im, re * rhs.im + im * rhs.re};
    }
    Complex16 operator*(Q15 scalar) const { return {re * scalar, im * scalar}; }
    Complex16 operator-() const { return {-re, -im}; }

    Complex16& operator+=(const Complex16& rhs) { return *this = *this + rhs; }
    Complex16& operator-=(const Complex16& rhs) { return *this = *this - rhs; }
    Complex16& operator*=(const Complex16& rhs) { return *this = *this * rhs; }

    // Complex conjugate.
    Complex16 conj() const { return {re, -im}; }
    // Magnitude squared: re^2 + im^2.
    Q15 magsq() const { return re * re + im * im; }

    constexpr bool operator==(const Complex16& rhs) const noexcept { return re == rhs.re && im == rhs.im; }
    constexpr bool operator!=(const Complex16& rhs) const noexcept { return !(*this == rhs); }
};

inline Complex16 operator*(Q15 scalar, const Complex16& c) { return c * scalar; }

} // namespace dsp

#endif
