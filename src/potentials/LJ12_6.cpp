//
// Created by alex on 8/4/24.
//

#include "LJ12_6.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"

LJ12_6::LJ12_6() : m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)) { }

void LJ12_6::handleCell(Cell &cell) {
    Kokkos::parallel_for("LJ12-6 - Cell", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {cell.getNumIndices(), cell.getNumIndices()}),
                         LJ12_6_Force(p_soa.f(), p_soa.id(), p_soa.r(), p_soa.sigma(), p_soa.epsilon(), cell.indices(), m_cutoff2));
}

void LJ12_6::handleCellPair(Cell &cell0, Cell &cell1, const math::d3& cell0_shift, const math::d3& cell1_shift) {
    Kokkos::parallel_for("LJ12-6 - Cell Pair", Kokkos::MDRangePolicy({0, 0}, {cell0.getNumIndices(), cell1.getNumIndices()}),
                         LJ12_6_ForcePair(p_soa.f(), p_soa.id(), p_soa.r(), p_soa.sigma(), p_soa.epsilon(), cell0.indices(), cell1.indices(), m_cutoff2, cell0_shift, cell1_shift));
}

void LJ12_6::LJ12_6_Force::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once
    const uint64_t s_idx_0 = indices[idx_0];
    const uint64_t s_idx_1 = indices[idx_1];
    
    if (id(s_idx_0) == id(s_idx_1)) return;

    const math::d3 dr = r[s_idx_0] - r[s_idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps[s_idx_0] * eps[s_idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const auto force = dr * fac;

    Kokkos::atomic_add(&f[s_idx_0], force);
    Kokkos::atomic_sub(&f[s_idx_1], force);
}

void LJ12_6::LJ12_6_ForcePair::operator()(int idx_0, int idx_1) const {
    const uint64_t s_idx_0 = indices0[idx_0];
    const uint64_t s_idx_1 = indices1[idx_1];
    
    if (id[s_idx_0] == id[s_idx_1]) return;

    const math::d3 dr = (r[s_idx_0] + shift0) - (r[s_idx_1] + shift1);
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps[s_idx_0] * eps[s_idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const auto force = dr * fac;

    Kokkos::atomic_add(&f[s_idx_0], force);
    Kokkos::atomic_sub(&f[s_idx_1], force);
}
