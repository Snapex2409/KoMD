//
// Created by alex on 8/4/24.
//

#include "Grav.h"
#include "molecule/Site.h"
#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"

Grav::Grav() { }

void Grav::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("Grav", pairList.size(), Grav_Force(pairList.getPairs(), pairList.getOffsets(), p_soa.f(), p_soa.id(), p_soa.r(), p_soa.mass()));
}

void Grav::Grav_Force::operator()(int idx) const {
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);
    
    if (id(s_idx_0) == id(s_idx_1)) return;
    // we ignore halo effects for gravity
    if (pair_offsets(idx, 0) != 0 || pair_offsets(idx, 1) != 0) return;

    const math::d3 dr = r[s_idx_0] - r[s_idx_1];
    const double dr2 = dr.dot(dr);
    const double fac = - m[s_idx_0] * m[s_idx_1] / (Kokkos::sqrt(dr2) * dr2);
    const auto force = dr * fac;

    Kokkos::atomic_add(&f[s_idx_0], force);
    Kokkos::atomic_sub(&f[s_idx_1], force);
}
