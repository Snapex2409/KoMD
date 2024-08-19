//
// Created by alex on 8/4/24.
//

#include "LJ12_6.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

LJ12_6::LJ12_6() : m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)) {}

void LJ12_6::handleCell(Cell &cell) {
    SOA& soa = cell.soa();
    Kokkos::parallel_for("LJ12-6 - Cell", Kokkos::MDRangePolicy({0, 0}, {soa.size(), soa.size()}), LJ12_6_Force(soa, m_cutoff2));
}

void LJ12_6::handleCellPair(Cell &cell0, Cell &cell1) {
    SOA& soa0 = cell0.soa();
    SOA& soa1 = cell1.soa();
    Kokkos::parallel_for("LJ12-6 - Cell Pair", Kokkos::MDRangePolicy({0, 0}, {soa0.size(), soa1.size()}), LJ12_6_ForcePair(soa0, soa1, m_cutoff2));
}

void LJ12_6::LJ12_6_Force::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once

    auto access = soa.fScatter().access();

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
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;

    access(idx_0) += dr * fac;
    access(idx_1) += dr * (-fac);
}

void LJ12_6::LJ12_6_ForcePair::operator()(int idx_0, int idx_1) const {
    auto access0 = soa0.fScatter().access();
    auto access1 = soa1.fScatter().access();

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
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;

    access0(idx_0) += dr * fac;
    access1(idx_1) += dr * (-fac);
}
