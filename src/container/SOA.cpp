//
// Created by alex on 7/31/24.
//

#include "SOA.h"

SOA::SOA() : m_current_size(0) { }

void SOA::resize(uint64_t newSize) {
    Kokkos::resize(m_epsilon, newSize);
    Kokkos::resize(m_sigma, newSize);
    Kokkos::resize(m_mass, newSize);
    Kokkos::resize(m_id, newSize);
    Kokkos::resize(m_idx, newSize);
    Kokkos::resize(m_r, newSize);
    Kokkos::resize(m_f, newSize);
    Kokkos::resize(m_v, newSize);

    m_current_size = newSize;
}

void SOA::createBuffers(uint64_t size) {
    m_epsilon = KW::vec_t<double>("SOA_eps", size);
    m_sigma = KW::vec_t<double>("SOA_sig", size);
    m_mass = KW::vec_t<double>("SOA_m", size);
    m_id = KW::vec_t<uint64_t>("SOA_id", size);
    m_idx = KW::vec_t<uint64_t>("SOA_idx", size);
    m_r = KW::vec_t<math::d3>("SOA_r", size);
    m_f = KW::vec_t<math::d3>("SOA_f", size);
    m_v = KW::vec_t<math::d3>("SOA_v", size);
    m_current_size = size;
}
