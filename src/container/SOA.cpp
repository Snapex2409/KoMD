//
// Created by alex on 7/31/24.
//

#include "SOA.h"

SOA::SOA(uint64_t size) : m_current_size(size) {
    m_epsilon.resize(size);
    m_sigma.resize(size);
    m_mass.resize(size);
    m_id.resize(size);
    m_r.resize(size);
    m_f.resize(size);
    m_fold.resize(size);
    m_v.resize(size);
}

void SOA::resize(uint64_t newSize) {
    m_epsilon.resize(newSize);
    m_sigma.resize(newSize);
    m_mass.resize(newSize);
    m_id.resize(newSize);
    m_r.resize(newSize);
    m_f.resize(newSize);
    m_fold.resize(newSize);
    m_v.resize(newSize);

    m_current_size = newSize;
}
