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
    Kokkos::resize(m_r, newSize);
    Kokkos::resize(m_f, newSize);
    Kokkos::resize(m_v, newSize);

    m_current_size = newSize;
}

void SOA::createBuffers() {
    m_epsilon = KW::vec_t<double>("SOA_eps", 1);
    m_sigma = KW::vec_t<double>("SOA_sig", 1);
    m_mass = KW::vec_t<double>("SOA_m", 1);
    m_id = KW::vec_t<uint64_t>("SOA_id", 1);
    m_r = KW::vec_t<math::d3>("SOA_r", 1);
    m_f = KW::vec_t<math::d3>("SOA_f", 1);
    m_v = KW::vec_t<math::d3>("SOA_v", 1);
}
