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
    SOA& soa = cell.soa();
    Kokkos::parallel_for("LJ12-6 - Cell", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {soa.size(), soa.size()}),
                         LJ12_6_Force(soa.f(), soa.id(), soa.r(), soa.sigma(), soa.epsilon(), m_cutoff2));
}

void LJ12_6::handleCellPair(Cell &cell0, Cell &cell1) {
    SOA& soa0 = cell0.soa();
    SOA& soa1 = cell1.soa();
    Kokkos::parallel_for("LJ12-6 - Cell Pair", Kokkos::MDRangePolicy({0, 0}, {soa0.size(), soa1.size()}),
                         LJ12_6_ForcePair(soa0.f(), soa1.f(), soa0.id(), soa1.id(), soa0.r(), soa1.r(), soa0.sigma(), soa1.sigma(), soa0.epsilon(), soa1.epsilon(), m_cutoff2));
}

void LJ12_6::LJ12_6_Force::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once
    if (id(idx_0) == id(idx_1)) return;

    const math::d3 dr = r[idx_0] - r[idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[idx_0] + sig[idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps[idx_0] * eps[idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const auto force = dr * fac;

    Kokkos::atomic_add(&f[idx_0], force);
    Kokkos::atomic_sub(&f[idx_1], force);
}

void LJ12_6::LJ12_6_ForcePair::operator()(int idx_0, int idx_1) const {
    if (id0[idx_0] == id1[idx_1]) return;

    const math::d3 dr = r0[idx_0] - r1[idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig0[idx_0] + sig1[idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps0[idx_0] * eps1[idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const auto force = dr * fac;

    Kokkos::atomic_add(&f0[idx_0], force);
    Kokkos::atomic_sub(&f1[idx_1], force);
}
