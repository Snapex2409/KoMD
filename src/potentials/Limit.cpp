//
// Created by alex on 8/4/24.
//

#include "Limit.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"
#include "math/Array.h"

Limit::Limit() : m_limit_factor(Registry::instance->configuration()->limit_factor) { }

void Limit::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("Limit - Cell", p_soa.size(),
                         Limit_Force(p_soa.f(), p_soa.sigma(), p_soa.epsilon(), m_limit_factor));
}

void Limit::Limit_Force::operator()(int idx) const {
    const uint64_t s_idx = idx;

    const double fmax = limit_factor * eps[s_idx] / sig[s_idx];
    const math::d3 force = f[s_idx];
    f[s_idx] = math::max(math::min(force, fmax), -fmax);
}
