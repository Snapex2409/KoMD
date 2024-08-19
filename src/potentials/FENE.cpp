//
// Created by alex on 8/4/24.
//

#include "FENE.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

FENE::FENE() : m_stiffness_factor(Registry::instance->configuration()->stiffness_factor) {
    p_run_pairs = false;
}

void FENE::handleCell(Cell &cell) {
    SOA& soa = cell.soa();
    Kokkos::parallel_for("FENE - Cell", Kokkos::MDRangePolicy({0, 0}, {soa.size(), soa.size()}), FENE_Force(soa, m_stiffness_factor));
}

void FENE::handleCellPair(Cell &cell0, Cell &cell1) { }

void FENE::FENE_Force::operator()(int idx_0, int idx_1) const {
    if (idx_1 <= idx_0) return; // only compute pair once

    auto access = soa.fScatter().access();

    auto& id = soa.id();
    auto& r = soa.r();
    auto& sig = soa.sigma();
    auto& eps = soa.epsilon();

    if (id[idx_0] != id[idx_1]) return;

    const math::d3 dr = r[idx_1] - r[idx_0];
    const math::d3 invdr = math::d3{1, 1, 1} / dr;
    const double dr2 = dr.dot(dr);

    const double sigma = (sig[idx_0] + sig[idx_1]) / 2.0;
    const double R02 = std::pow(1.5 * sigma, 2);
    if (R02 <= dr2) {
        access(idx_0) += invdr * 1e+30;
        access(idx_1) += invdr * (-1e+30);
        return;
    }

    const double epsilon = std::sqrt(eps[idx_0] * eps[idx_1]);
    const double k = stiffness_factor * epsilon / std::pow(sigma, 2);

    const double fac = k / (1 - (dr2/R02));

    access(idx_0) += dr * fac;
    access(idx_1) += dr * (-fac);
}
