//
// Created by alex on 8/4/24.
//

#ifndef KOMD_SCIVALUE_H
#define KOMD_SCIVALUE_H

#include "Kokkos_Core.hpp"

/**
 * Represents a floating point number in scientific notation
 * */
class SciValue {
public:
    constexpr SciValue(double factor, int64_t exp) : m_factor(factor), m_exp(exp) {
        normalize();
    }
    constexpr explicit SciValue(double value) : m_factor(value), m_exp(0) {
        normalize();
    }

    KOKKOS_FUNCTION SciValue& operator*=(const SciValue& other);
    KOKKOS_FUNCTION SciValue& operator/=(const SciValue& other);
    KOKKOS_FUNCTION SciValue& operator*=(double other);
    KOKKOS_FUNCTION SciValue& operator/=(double other);

    KOKKOS_FUNCTION SciValue operator*(const SciValue& other) const;
    KOKKOS_FUNCTION SciValue operator/(const SciValue& other) const;
    KOKKOS_FUNCTION SciValue operator*(double other) const;
    KOKKOS_FUNCTION SciValue operator/(double other) const;

    KOKKOS_FUNCTION operator double() const; // NOLINT(*-explicit-constructor)
private:
    /**
     * normalizes this number such that the absolute value of the mantissa is < 10
     * */
    KOKKOS_FUNCTION constexpr void normalize() {
        double abs_val = std::abs(m_factor);
        while (abs_val >= 10.0) {
            m_factor /= 10.0;
            m_exp += 1;
            abs_val = std::abs(m_factor);
        }
        while (abs_val <= 0.1) {
            m_factor *= 10.0;
            m_exp -= 1;
            abs_val = std::abs(m_factor);
        }
    }
    /// "mantissa"
    double m_factor;
    /// "exponent"
    int64_t m_exp;
};


#endif //KOMD_SCIVALUE_H
