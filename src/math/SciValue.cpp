//
// Created by alex on 8/4/24.
//

#include <cmath>
#include "SciValue.h"

SciValue::operator double() const {
    return m_factor * (double) m_exp;
}

SciValue &SciValue::operator*=(const SciValue &other) {
    m_factor *= other.m_factor;
    m_exp += other.m_exp;

    normalize();
    return *this;
}

SciValue &SciValue::operator/=(const SciValue &other) {
    m_factor /= other.m_factor;
    m_exp -= other.m_exp;

    normalize();
    return *this;
}

SciValue &SciValue::operator*=(double other) {
    SciValue sv (other);
    return operator*=(sv);
}

SciValue &SciValue::operator/=(double other) {
    SciValue sv (other);
    return operator/=(sv);
}

SciValue SciValue::operator*(const SciValue &other) const {
    return {m_factor*other.m_factor, m_exp+other.m_exp};
}

SciValue SciValue::operator/(const SciValue &other) const {
    return {m_factor/other.m_factor, m_exp-other.m_exp};
}

SciValue SciValue::operator*(double other) const {
    return operator*(SciValue(other));
}

SciValue SciValue::operator/(double other) const {
    return operator/(SciValue(other));
}
