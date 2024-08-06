//
// Created by alex on 8/4/24.
//

#include "Limit.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"

Limit::Limit() : m_limit_factor(Registry::instance->configuration()->limit_factor) {
    p_run_pairs = false;
}

void Limit::handleCell(Cell &cell) {
    if (!p_use_soa) {
        auto& molecules = cell.molecules();
        for (uint64_t mi = 0; mi < molecules.size(); mi++) {
            Molecule& mol_i = molecules[mi];
            uint64_t num_sites = mol_i.getSites().size();
            for (uint64_t si = 0; si < num_sites; si++) {
                Site& site = mol_i.getSites()[si];
                computeForce(site);
            }
        }
    }
        // use SOA
    else {
        auto& soa = cell.soa();
        auto size = soa.size();
        for (uint64_t idx_i = 0; idx_i < size; idx_i++) {
            computeForceSOA(idx_i, soa.r(), soa.f(), soa.sigma(), soa.epsilon());
        }
    }
}

void Limit::handleCellPair(Cell &cell0, Cell &cell1) { }

void Limit::computeForce(Site &site) const {
    const double fmax = m_limit_factor * site.getEpsilon() / site.getSigma();
    const math::d3 f = site.f_arr();
    site.f_arr() = math::max(math::min(f, fmax), -fmax);
}

void Limit::computeForceSOA(uint64_t idx, SOA::vec_t<math::d3>& r0, SOA::vec_t<math::d3>& f0, SOA::vec_t<double>& sigmas, SOA::vec_t<double>& epsilons) const {
    const double fmax = m_limit_factor * epsilons[idx] / sigmas[idx];
    const math::d3 f = f0[idx];
    f0[idx] = math::max(math::min(f, fmax), -fmax);
}
