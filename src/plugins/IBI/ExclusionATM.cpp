//
// Created by alex on 11/11/24.
//

#include "ExclusionATM.h"

#include <math/Geometry.h>

#include "Registry.h"

ExclusionATM::ExclusionATM() :
    m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)),
    m_nu(Registry::instance->configuration()->energy_3b),
    m_exclusion_low(Registry::instance->configuration()->IBI_exclusion_low),
    m_exclusion_high(Registry::instance->configuration()->IBI_exclusion_high) { }

void ExclusionATM::handleTripleList(TripleList &tripleList) {
    Kokkos::parallel_for("ATM", tripleList.size(),
        ExATM_Force(tripleList.getTriplets(), tripleList.getOffsets(),
            p_soa.f(), p_soa.r(), m_nu, m_cutoff2, m_exclusion_low, m_exclusion_high));
}

void ExclusionATM::ExATM_Force::operator()(int idx) const {
    const uint64_t s_idx_0 = triplets(idx, 0);
    const uint64_t s_idx_1 = triplets(idx, 1);
    const uint64_t s_idx_2 = triplets(idx, 2);

    const math::d3 shift0 = triple_offsets(idx, 0);
    const math::d3 shift1 = triple_offsets(idx, 1);
    const math::d3 shift2 = triple_offsets(idx, 2);

    const math::d3 r0 = r[s_idx_0] + shift0;
    const math::d3 r1 = r[s_idx_1] + shift1;
    const math::d3 r2 = r[s_idx_2] + shift2;

    if (!math::pointInBox(r0, exclusion_low, exclusion_high) &&
        !math::pointInBox(r1, exclusion_low, exclusion_high) &&
        !math::pointInBox(r2, exclusion_low, exclusion_high)) return;

    const math::d3 dr0_1 = r0 - r1;
    const math::d3 dr0_2 = r0 - r2;
    const math::d3 dr1_2 = r1 - r2;

    const double dr0_1_2 = dr0_1.dot(dr0_1);
    const double dr0_2_2 = dr0_2.dot(dr0_2);
    const double dr1_2_2 = dr1_2.dot(dr1_2);
    // symmetrical condition
    if (dr0_1_2 > cutoff2 || dr0_2_2 > cutoff2 || dr1_2_2 > cutoff2) return;

    // we are within range
    // compute forces
    const double r0_1 = Kokkos::sqrt(dr0_1_2);
    const double r0_2 = Kokkos::sqrt(dr0_2_2);
    const double r1_2 = Kokkos::sqrt(dr1_2_2);

    const double dVdR01 = comp_force(nu, r0_1, r0_2, r1_2);
    const double dVdR02 = comp_force(nu, r0_2, r0_1, r1_2);

    const math::d3 f0 = dr0_1 * dVdR01 + dr0_2 * dVdR02;
    const math::d3 f1 = dr0_1 * (-dVdR01) + dr0_2 * dVdR02;
    const math::d3 f2 = dr0_1 * (-dVdR01) + dr0_2 * (-dVdR02);

    if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], f0);
    if (shift1 == 0) Kokkos::atomic_add(&f[s_idx_1], f1);
    if (shift2 == 0) Kokkos::atomic_add(&f[s_idx_2], f2);
}