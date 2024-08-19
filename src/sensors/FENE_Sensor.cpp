//
// Created by alex on 8/8/24.
//

#include "FENE_Sensor.h"
#include "Registry.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

FENE_Sensor::FENE_Sensor() :
Potential_Sensor("FENE", Registry::instance->configuration()->sensor_fene_bins),
m_stiffness_factor(Registry::instance->configuration()->stiffness_factor) {
    p_run_pairs = false;
}

void FENE_Sensor::handleCell(Cell &cell) {
    SOA& soa = cell.soa();
    Kokkos::parallel_for("Sensor - FENE - Cell", Kokkos::MDRangePolicy({0, 0}, {soa.size(), soa.size()}), FENE_Pot(soa, *this, m_stiffness_factor, p_sigma, p_bins));
}

void FENE_Sensor::handleCellPair(Cell &cell0, Cell &cell1) { }

void FENE_Sensor::FENE_Pot::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once

    auto data_f_access = sensor.p_data_f_scatter.access();
    auto data_u_access = sensor.p_data_u_scatter.access();
    auto count_f_access = sensor.p_count_f_scatter.access();
    auto count_u_access = sensor.p_count_u_scatter.access();


    auto& id = soa.id();
    auto& r = soa.r();
    auto& sig = soa.sigma();
    auto& eps = soa.epsilon();

    if (id[idx_0] != id[idx_1]) return;

    const math::d3 dr = r[idx_1] - r[idx_0];
    const double dr2 = dr.dot(dr);

    const double sigma = (sig[idx_0] + sig[idx_1]) / 2.0;
    const double R02 = std::pow(1.5 * sigma, 2);
    if (R02 <= dr2) return;

    const double epsilon = std::sqrt(eps[idx_0] * eps[idx_1]);
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
