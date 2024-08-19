//
// Created by alex on 8/8/24.
//

#include "LJ12_6_Sensor.h"
#include "Registry.h"

LJ12_6_Sensor::LJ12_6_Sensor() :
Potential_Sensor("LJ12_6", Registry::instance->configuration()->sensor_lj_bins),
m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)), m_total_pot("Sensor Total Pot", 1) {
    m_total_pot_scatter = Kokkos::Experimental::ScatterView<double*>(m_total_pot);
}

void LJ12_6_Sensor::measure() {
    m_total_pot[0] = 0;
    Potential_Sensor::measure();
    Kokkos::Experimental::contribute(m_total_pot, m_total_pot_scatter);
}

void LJ12_6_Sensor::handleCell(Cell &cell) {
    SOA& soa = cell.soa();
    Kokkos::parallel_for("Sensor - LJ12-6 - Cell", Kokkos::MDRangePolicy({0, 0}, {soa.size(), soa.size()}), LJ12_6_Pot(soa, *this, m_cutoff2, p_sigma, p_bins));
}

void LJ12_6_Sensor::handleCellPair(Cell &cell0, Cell &cell1) {
    SOA& soa0 = cell0.soa();
    SOA& soa1 = cell1.soa();
    Kokkos::parallel_for("Sensor - LJ12-6 - Cell Pair", Kokkos::MDRangePolicy({0, 0}, {soa0.size(), soa1.size()}), LJ12_6_PotPair(soa0, soa1, *this, m_cutoff2, p_sigma, p_bins));
}

void LJ12_6_Sensor::LJ12_6_Pot::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once

    auto data_f_access = sensor.p_data_f_scatter.access();
    auto data_u_access = sensor.p_data_u_scatter.access();
    auto count_f_access = sensor.p_count_f_scatter.access();
    auto count_u_access = sensor.p_count_u_scatter.access();
    auto total_pot_access = sensor.m_total_pot_scatter.access();

    auto& id = soa.id();
    auto& r = soa.r();
    auto& sig = soa.sigma();
    auto& eps = soa.epsilon();

    if (id[idx_0] == id[idx_1]) return;

    const math::d3 dr = r[idx_0] - r[idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[idx_0] + sig[idx_1]) / 2.0;
    const double epsilon = std::sqrt(eps[idx_0] * eps[idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = std::pow(sig2 * invdr2, 3.0);
    const double lj12 = std::pow(lj6, 2.0);
    const double r_pos = dr.L2();

    const double f = 24.0 * epsilon * (2.0 * lj12 - lj6) / r_pos;
    const double u = 4.0 * epsilon * (lj12 - lj6);

    auto bin = get_bin(r_pos, max_sigma, bins);
    data_f_access(bin) += f;
    data_u_access(bin) += u;
    count_f_access(bin) += 1;
    count_u_access(bin) += 1;
    total_pot_access(0) += u;
}

void LJ12_6_Sensor::LJ12_6_PotPair::operator()(int idx_0, int idx_1) const {
    auto data_f_access = sensor.p_data_f_scatter.access();
    auto data_u_access = sensor.p_data_u_scatter.access();
    auto count_f_access = sensor.p_count_f_scatter.access();
    auto count_u_access = sensor.p_count_u_scatter.access();
    auto total_pot_access = sensor.m_total_pot_scatter.access();

    auto& id0 = soa0.id();
    auto& id1 = soa1.id();
    auto& r0 = soa0.r();
    auto& r1 = soa1.r();
    auto& sig0 = soa0.sigma();
    auto& sig1 = soa1.sigma();
    auto& eps0 = soa0.epsilon();
    auto& eps1 = soa1.epsilon();

    if (id0[idx_0] == id1[idx_1]) return;

    const math::d3 dr = r0[idx_0] - r1[idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig0[idx_0] + sig1[idx_1]) / 2.0;
    const double epsilon = std::sqrt(eps0[idx_0] * eps1[idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = std::pow(sig2 * invdr2, 3.0);
    const double lj12 = std::pow(lj6, 2.0);
    const double r_pos = dr.L2();

    const double f = 24.0 * epsilon * (2.0 * lj12 - lj6) / r_pos;
    const double u = 4.0 * epsilon * (lj12 - lj6);

    auto bin = get_bin(r_pos, max_sigma, bins);
    data_f_access(bin) += f;
    data_u_access(bin) += u;
    count_f_access(bin) += 1;
    count_u_access(bin) += 1;
    total_pot_access(0) += u;
}
