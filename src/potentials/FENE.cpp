//
// Created by alex on 8/4/24.
//

#include "FENE.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"

FENE::FENE() : m_stiffness_factor(Registry::instance->configuration()->stiffness_factor) {
    p_run_pairs = false;
}

void FENE::handleCell(Cell &cell) {
    Kokkos::parallel_for("FENE - Cell", Kokkos::MDRangePolicy({0, 0}, {cell.getNumIndices(), cell.getNumIndices()}), 
                         FENE_Force(p_soa.f(), p_soa.id(), p_soa.r(), p_soa.sigma(), p_soa.epsilon(), cell.indices(), m_stiffness_factor));
}

void FENE::handleCellPair(Cell &cell0, Cell &cell1, const math::d3& cell1_shift) { }

void FENE::FENE_Force::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once
    const uint64_t s_idx_0 = indices[idx_0];
    const uint64_t s_idx_1 = indices[idx_1];
    
    if (id[s_idx_0] != id[s_idx_1]) return;

    const math::d3 dr = r[s_idx_1] - r[s_idx_0];
    const math::d3 invdr = math::d3{1, 1, 1} / dr;
    const double dr2 = dr.dot(dr);

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double R02 = std::pow(1.5 * sigma, 2);
    if (R02 <= dr2) {
        const auto force = invdr * 1e+30;
        Kokkos::atomic_add(&f[s_idx_0], force);
        Kokkos::atomic_sub(&f[s_idx_1], force);
        return;
    }

    const double epsilon = std::sqrt(eps[s_idx_0] * eps[s_idx_1]);
    const double k = stiffness_factor * epsilon / std::pow(sigma, 2);

    const double fac = k / (1 - (dr2/R02));
    const auto force = dr * fac;

    Kokkos::atomic_add(&f[s_idx_0], force);
    Kokkos::atomic_sub(&f[s_idx_1], force);
}
