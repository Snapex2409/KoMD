//
// Created by alex on 8/4/24.
//

#include "Limit.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"
#include "math/Array.h"

#include "Kokkos_Core.hpp"

Limit::Limit() : m_limit_factor(Registry::instance->configuration()->limit_factor) {
    p_run_pairs = false;
    p_run_contribution = false;
}

void Limit::handleCell(Cell &cell) {
    SOA& soa = cell.soa();
    Kokkos::parallel_for("Limit - Cell", soa.size(), Limit_Force(soa, m_limit_factor));
}

void Limit::handleCellPair(Cell &cell0, Cell &cell1) { }

void Limit::Limit_Force::operator()(int idx) const {
    auto& eps = soa.epsilon();
    auto& sig = soa.sigma();
    auto& f = soa.f();

    const double fmax = limit_factor * eps[idx] / sig[idx];
    const math::d3 force = f[idx];
    f[idx] = math::max(math::min(force, fmax), -fmax);
}
