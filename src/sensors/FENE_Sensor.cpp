//
// Created by alex on 8/8/24.
//

#include "FENE_Sensor.h"
#include "Registry.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

FENE_Sensor::FENE_Sensor() :
Potential_Sensor("FENE", Registry::instance->configuration()->sensor_fene_bins),
m_stiffness_factor(Registry::instance->configuration()->stiffness_factor) { }

void FENE_Sensor::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("FENE", pairList.size(),
                         FENE_Pot(pairList.getPairs(), pairList.getOffsets(), p_soa.id(), p_soa.r(), p_soa.sigma(), p_soa.epsilon(), p_data_u_scatter, p_data_f_scatter, p_count_u_scatter, p_count_f_scatter,
                                  m_stiffness_factor, p_max_sigma, p_bins));
}

void FENE_Sensor::FENE_Pot::operator()(int idx) const {
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);
    
    if (id[s_idx_0] != id[s_idx_1]) return;
    if (pair_offsets(idx, 0) != 0 || pair_offsets(idx, 1) != 0) return;

    auto data_f_access = data_f_scatter.access();
    auto data_u_access = data_u_scatter.access();
    auto count_f_access = count_f_scatter.access();
    auto count_u_access = count_u_scatter.access();

    const math::d3 dr = r[s_idx_1] - r[s_idx_0];
    const double dr2 = dr.dot(dr);

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double R02 = std::pow(1.5 * sigma, 2);
    if (R02 <= dr2) return;

    const double epsilon = std::sqrt(eps[s_idx_0] * eps[s_idx_1]);
    const double k = stiffness_factor * epsilon / std::pow(sigma, 2);
    const double r_pos = dr.L2();

    const double f = r_pos * k / (1 - (dr2/R02));
    const double u = -0.5 * k * R02 * std::log(1 - (dr2/R02));

    auto bin = get_bin(r_pos, max_sigma, bins);
    data_f_access(bin) += f;
    data_u_access(bin) += u;
    count_f_access(bin) += 1;
    count_u_access(bin) += 1;
}
