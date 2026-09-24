#include "../dsptypes.hpp"

#include <cmath>

namespace dsp {

Q15 Q15::operator+(Q15 rhs) const {
    int32_t sum = static_cast<int32_t>(raw_) + rhs.raw_;
    if (sum > MAX) sum = MAX;
    if (sum < MIN) sum = MIN;
    return Q15(static_cast<int16_t>(sum));
}

Q15 Q15::operator-(Q15 rhs) const {
    int32_t diff = static_cast<int32_t>(raw_) - rhs.raw_;
    if (diff > MAX) diff = MAX;
    if (diff < MIN) diff = MIN;
    return Q15(static_cast<int16_t>(diff));
}

Q15 Q15::operator*(Q15 rhs) const {
    int32_t product = static_cast<int32_t>(raw_) * rhs.raw_;
    int32_t rounded = (product + (1 << 14)) >> 15;
    if (rounded > MAX) rounded = MAX;
    if (rounded < MIN) rounded = MIN;
    return Q15(static_cast<int16_t>(rounded));
}

Q15 Q15::operator/(Q15 rhs) const {
    if (rhs.raw_ == 0) {
        return raw_ >= 0 ? Q15(MAX) : Q15(MIN);
    }
    int64_t numerator = static_cast<int64_t>(raw_) << 15;
    int64_t quotient = numerator / rhs.raw_;
    if (quotient > MAX) quotient = MAX;
    if (quotient < MIN) quotient = MIN;
    return Q15(static_cast<int16_t>(quotient));
}

Q15 Q15::operator-() const {
    return Q15() - *this;
}

Q15 Q15::from_double(double x) {
    double scaled = x * 32768.0;
    if (scaled >= MAX) return Q15(MAX);
    if (scaled <= MIN) return Q15(MIN);
    return Q15(static_cast<int16_t>(std::lround(scaled)));
}

Q15 Q15::from_float(float x) {
    return from_double(static_cast<double>(x));
}

double Q15::to_double() const {
    return static_cast<double>(raw_) / 32768.0;
}

float Q15::to_float() const {
    return static_cast<float>(raw_) / 32768.0f;
}

} // namespace dsp
