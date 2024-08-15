//
// Created by alex on 7/31/24.
//

#include "SOA.h"

SOA::SOA(uint64_t size) : m_current_size(size) {
    m_epsilon = vec_t<double>("SOA_eps", size);
    m_sigma = vec_t<double>("SOA_", size);
    m_mass = vec_t<double>("SOA_", size);
    m_id = vec_t<uint64_t>("SOA_", size);
    m_r = vec_t<math::d3>("SOA_", size);
    m_f = vec_t<math::d3>("SOA_", size);
    m_v = vec_t<math::d3>("SOA_", size);
}

void SOA::resize(uint64_t newSize) {
    m_epsilon = vec_t<double>("SOA_eps", newSize);
    m_sigma = vec_t<double>("SOA_", newSize);
    m_mass = vec_t<double>("SOA_", newSize);
    m_id = vec_t<uint64_t>("SOA_", newSize);
    m_r = vec_t<math::d3>("SOA_", newSize);
    m_f = vec_t<math::d3>("SOA_", newSize);
    m_v = vec_t<math::d3>("SOA_", newSize);

    m_current_size = newSize;
}
