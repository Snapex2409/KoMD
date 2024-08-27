//
// Created by alex on 8/8/24.
//

#include "LJ12_6_Sensor.h"
#include "Registry.h"

LJ12_6_Sensor::LJ12_6_Sensor() :
Potential_Sensor("LJ12_6", Registry::instance->configuration()->sensor_lj_bins),
m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)), m_total_pot("Sensor Total Pot", 1) {
    m_total_pot_scatter = KW::vec_scatter_t<double>(m_total_pot);
}

void LJ12_6_Sensor::measure() {
    m_total_pot[0] = 0;
    m_total_pot_scatter.reset();
    Potential_Sensor::measure();
    Kokkos::Experimental::contribute(m_total_pot, m_total_pot_scatter);
}

void LJ12_6_Sensor::handleCell(Cell &cell) {
    Kokkos::parallel_for("Sensor - LJ12-6 - Cell", Kokkos::MDRangePolicy({0, 0}, {cell.getNumIndices(), cell.getNumIndices()}),
                         LJ12_6_Pot(p_soa.id(), p_soa.r(), p_soa.sigma(), p_soa.epsilon(), cell.indices(),
                                    p_data_u_scatter, p_data_f_scatter, p_count_u_scatter, p_count_f_scatter,
                                    m_total_pot_scatter, m_cutoff2, p_max_sigma, p_bins));
}

void LJ12_6_Sensor::handleCellPair(Cell &cell0, Cell &cell1, const math::d3& cell1_shift) {
    Kokkos::parallel_for("Sensor - LJ12-6 - Cell Pair", Kokkos::MDRangePolicy({0, 0}, {cell0.getNumIndices(), cell1.getNumIndices()}),
                         LJ12_6_PotPair(p_soa.id(), p_soa.r(), p_soa.sigma(), p_soa.epsilon(), cell0.indices(), cell1.indices(),
                                        p_data_u_scatter, p_data_f_scatter, p_count_u_scatter, p_count_f_scatter,
                                        m_total_pot_scatter, m_cutoff2, p_max_sigma, p_bins, cell1_shift));
}

void LJ12_6_Sensor::LJ12_6_Pot::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once
    const uint64_t s_idx_0 = indices[idx_0];
    const uint64_t s_idx_1 = indices[idx_1];
    
    if (id[s_idx_0] == id[s_idx_1]) return;

    auto data_f_access = data_f_scatter.access();
    auto data_u_access = data_u_scatter.access();
    auto count_f_access = count_f_scatter.access();
    auto count_u_access = count_u_scatter.access();
    auto total_pot_access = total_pot_scatter.access();

    const math::d3 dr = r[s_idx_0] - r[s_idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double epsilon = std::sqrt(eps[s_idx_0] * eps[s_idx_1]);

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
    auto data_f_access = data_f_scatter.access();
    auto data_u_access = data_u_scatter.access();
    auto count_f_access = count_f_scatter.access();
    auto count_u_access = count_u_scatter.access();
    auto total_pot_access = total_pot_scatter.access();

    const uint64_t s_idx_0 = indices0[idx_0];
    const uint64_t s_idx_1 = indices1[idx_1];
    if (id[s_idx_0] == id[s_idx_1]) return;

    const math::d3 dr = r[s_idx_0] - (r[s_idx_1] + shift);
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double epsilon = std::sqrt(eps[s_idx_0] * eps[s_idx_1]);

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
