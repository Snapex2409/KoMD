//
// Created by alex on 8/8/24.
//

#include "LJ12_6_Sensor.h"
#include "Registry.h"

LJ12_6_Sensor::LJ12_6_Sensor() :
Potential_Sensor("LJ12_6", Registry::instance->configuration()->sensor_lj_bins),
m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)) {}

void LJ12_6_Sensor::handleCell(Cell &cell) {
    if (!p_use_soa) {
        auto& molecules = cell.molecules();
        for (uint64_t mi = 0; mi < molecules.size(); mi++) {
            for (uint64_t mj = mi+1; mj < molecules.size(); mj++) {
                Molecule& mol_i = molecules[mi];
                Molecule& mol_j = molecules[mj];
                for (Site& si : mol_i.getSites()) {
                    for (Site& sj : mol_j.getSites()) {
                        computeForce(si, sj);
                    }
                }
            }
        }
    }
        // use SOA
    else {
        auto& soa = cell.soa();
        auto size = soa.size();
        for (uint64_t idx_i = 0; idx_i < size; idx_i++) {
            for (uint64_t idx_j = idx_i+1; idx_j < size; idx_j++) {
                if (soa.id()[idx_i] == soa.id()[idx_j]) continue;
                computeForceSOA(idx_i, idx_j, soa.r(), soa.r(), soa.f(), soa.f(), soa.sigma(), soa.sigma(), soa.epsilon(), soa.epsilon());
            }
        }
    }
}

void LJ12_6_Sensor::handleCellPair(Cell &cell0, Cell &cell1) {
    if (!p_use_soa) {
        auto& molecules0 = cell0.molecules();
        auto& molecules1 = cell1.molecules();
        for (uint64_t mi = 0; mi < molecules0.size(); mi++) {
            for (uint64_t mj = 0; mj < molecules1.size(); mj++) {
                Molecule& mol_i = molecules0[mi];
                Molecule& mol_j = molecules1[mj];
                for (Site& si : mol_i.getSites()) {
                    for (Site& sj : mol_j.getSites()) {
                        computeForce(si, sj);
                    }
                }
            }
        }
    }
        // use SOA
    else {
        auto& soa0 = cell0.soa();
        auto& soa1 = cell1.soa();
        auto size0 = soa0.size();
        auto size1 = soa1.size();
        for (uint64_t idx_i = 0; idx_i < size0; idx_i++) {
            for (uint64_t idx_j = 0; idx_j < size1; idx_j++) {
                computeForceSOA(idx_i, idx_j, soa0.r(), soa1.r(), soa0.f(), soa1.f(), soa0.sigma(), soa1.sigma(), soa0.epsilon(), soa1.epsilon());
            }
        }
    }
}

void LJ12_6_Sensor::computeForce(Site &site0, Site &site1) {
    const math::d3 dr = site0.r_arr() - site1.r_arr();
    const double dr2 = dr.dot(dr);
    if (dr2 > m_cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (site0.getSigma() + site1.getSigma()) / 2.0;
    const double epsilon = std::sqrt(site0.getEpsilon() * site1.getEpsilon());

    const double sig2 = sigma * sigma;
    const double lj6 = std::pow(sig2 * invdr2, 3.0);
    const double lj12 = std::pow(lj6, 2.0);
    const double r = dr.L2();

    const double f = 24.0 * epsilon * (2.0 * lj12 - lj6) / r;
    const double u = 4.0 * epsilon * (lj12 - lj6);

    auto bin = get_bin(r);
    p_data_f[bin] += f;
    p_data_u[bin] += u;
    p_count_f[bin] += 1;
    p_count_u[bin] += 1;
}

void LJ12_6_Sensor::computeForceSOA(uint64_t idx_0, uint64_t idx_1, SOA::vec_t<math::d3> &r0, SOA::vec_t<math::d3> &r1,
                                    SOA::vec_t<math::d3> &f0, SOA::vec_t<math::d3> &f1, SOA::vec_t<double> &sigmas0,
                                    SOA::vec_t<double> &sigmas1, SOA::vec_t<double> &epsilons0,
                                    SOA::vec_t<double> &epsilons1) {
    const math::d3 dr = r0[idx_0] - r1[idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > m_cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sigmas0[idx_0] + sigmas1[idx_1]) / 2.0;
    const double epsilon = std::sqrt(epsilons0[idx_0] * epsilons1[idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = std::pow(sig2 * invdr2, 3.0);
    const double lj12 = std::pow(lj6, 2.0);
    const double r = dr.L2();

    const double f = 24.0 * epsilon * (2.0 * lj12 - lj6) / r;
    const double u = 4.0 * epsilon * (lj12 - lj6);

    auto bin = get_bin(r);
    p_data_f[bin] += f;
    p_data_u[bin] += u;
    p_count_f[bin] += 1;
    p_count_u[bin] += 1;
}
