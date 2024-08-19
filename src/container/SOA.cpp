//
// Created by alex on 7/31/24.
//

#include "SOA.h"

SOA::SOA(uint64_t size) { resize(size); }

void SOA::resize(uint64_t newSize) {
    m_epsilon = vec_t<double>("SOA_eps", newSize);
    m_sigma = vec_t<double>("SOA_sig", newSize);
    m_mass = vec_t<double>("SOA_m", newSize);
    m_id = vec_t<uint64_t>("SOA_id", newSize);
    m_r = vec_t<math::d3>("SOA_r", newSize);
    m_f = vec_t<math::d3>("SOA_f", newSize);
    m_v = vec_t<math::d3>("SOA_v", newSize);

    m_epsilon_scatter = vec_scatter_t<double>(m_epsilon);
    m_sigma_scatter = vec_scatter_t<double>(m_sigma);
    m_mass_scatter = vec_scatter_t<double>(m_mass);
    m_id_scatter = vec_scatter_t<uint64_t>(m_id);
    m_r_scatter = vec_scatter_t<math::d3>(m_r);
    m_f_scatter = vec_scatter_t<math::d3>(m_f);
    m_v_scatter = vec_scatter_t<math::d3>(m_v);

    m_current_size = newSize;
}
