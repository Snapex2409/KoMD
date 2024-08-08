//
// Created by alex on 8/8/24.
//

#include "FENE_Sensor.h"
#include "Registry.h"

FENE_Sensor::FENE_Sensor() :
Potential_Sensor("FENE", Registry::instance->configuration()->sensor_fene_bins),
m_stiffness_factor(Registry::instance->configuration()->stiffness_factor) {
    p_run_pairs = false;
}

void FENE_Sensor::handleCell(Cell &cell) {
    if (!p_use_soa) {
        auto& molecules = cell.molecules();
        for (uint64_t mi = 0; mi < molecules.size(); mi++) {
            Molecule& mol_i = molecules[mi];
            uint64_t num_sites = mol_i.getSites().size();
            for (uint64_t si = 0; si < num_sites; si++) {
                for (uint64_t sj = si+1; sj < num_sites; sj++) {
                    Site& site_i = mol_i.getSites()[si];
                    Site& site_j = mol_i.getSites()[sj];
                    computeForce(site_i, site_j);
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
                if (soa.id()[idx_i] != soa.id()[idx_j]) break;
                computeForceSOA(idx_i, idx_j, soa.r(), soa.r(), soa.f(), soa.f(), soa.sigma(), soa.sigma(), soa.epsilon(), soa.epsilon());
            }
        }
    }
}

void FENE_Sensor::handleCellPair(Cell &cell0, Cell &cell1) { }

void FENE_Sensor::computeForce(Site &site0, Site &site1) {
    const math::d3 dr = site1.r_arr() - site0.r_arr();
    const double dr2 = dr.dot(dr);

    const double sigma = (site0.getSigma() + site1.getSigma()) / 2.0;
    const double R02 = std::pow(1.5 * sigma, 2);
    if (R02 <= dr2) return;

    const double epsilon = std::sqrt(site0.getEpsilon() * site1.getEpsilon());
    const double k = m_stiffness_factor * epsilon / std::pow(sigma, 2);
    const double r = dr.L2();

    const double f = r * k / (1 - (dr2/R02));
    const double u = -0.5 * k * R02 * std::log(1 - (dr2/R02));

    auto bin = get_bin(r);
    p_data_f[bin] += f;
    p_data_u[bin] += u;
    p_count_f[bin] += 1;
    p_count_u[bin] += 1;
}

void FENE_Sensor::computeForceSOA(uint64_t idx_0, uint64_t idx_1, SOA::vec_t<math::d3> &r0, SOA::vec_t<math::d3> &r1,
                                  SOA::vec_t<math::d3> &f0, SOA::vec_t<math::d3> &f1, SOA::vec_t<double> &sigmas0,
                                  SOA::vec_t<double> &sigmas1, SOA::vec_t<double> &epsilons0,
                                  SOA::vec_t<double> &epsilons1) {
    const math::d3 dr = r1[idx_1] - r0[idx_0];
    const double dr2 = dr.dot(dr);

    const double sigma = (sigmas0[idx_0] + sigmas1[idx_1]) / 2.0;
    const double R02 = std::pow(1.5 * sigma, 2);
    if (R02 <= dr2) return;

    const double epsilon = std::sqrt(epsilons0[idx_0] * epsilons1[idx_1]);
    const double k = m_stiffness_factor * epsilon / std::pow(sigma, 2);
    const double r = dr.L2();

    const double f = r * k / (1 - (dr2/R02));
    const double u = -0.5 * k * R02 * std::log(1 - (dr2/R02));

    auto bin = get_bin(r);
    p_data_f[bin] += f;
    p_data_u[bin] += u;
    p_count_f[bin] += 1;
    p_count_u[bin] += 1;
}
