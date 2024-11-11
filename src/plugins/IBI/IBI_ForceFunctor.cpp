//
// Created by alex on 11/5/24.
//

#include "IBI_ForceFunctor.h"
#include "Registry.h"
#include "math/Geometry.h"

//========================================
// Default IBI implementation
//========================================
template<>
IBI_ForceFunctor<IBI_Default>::IBI_ForceFunctor() : m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)),
    m_force(Registry::instance->configuration()->IBI_bins),
    m_potential(Registry::instance->configuration()->IBI_bins) { }

template<>
void IBI_ForceFunctor<IBI_Default>::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("IBI force", pairList.size(),
        IBI_Force(pairList.getPairs(), pairList.getOffsets(), p_soa.f(), p_soa.r(), m_force.getXValues(), m_force.getYValues(), m_force.getLowerDefault(), m_force.getUpperDefault(), m_cutoff2, m_exclusion_low, m_exclusion_high));
}

template<>
void IBI_ForceFunctor<IBI_Default>::IBI_Force::operator()(int idx) const {
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);
    const math::d3 shift0 = pair_offsets(idx, 0);
    const math::d3 shift1 = pair_offsets(idx, 1);

    const math::d3 r0 = r[s_idx_0] + shift0;
    const math::d3 r1 = r[s_idx_1] + shift1;
    const math::d3 dr = r0 - r1;
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;

    const auto r = Kokkos::sqrt(dr2);
    const auto F = FunctionPL::evaluateAt(r, force_def_low, force_def_high, force_x_values, force_y_values);
    const math::d3 force = (dr / r) * F;

    if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], force);
    if (shift1 == 0) Kokkos::atomic_sub(&f[s_idx_1], force);
}

//========================================
// Reload IBI implementation
//========================================
template<>
IBI_ForceFunctor<IBI_Reload>::IBI_ForceFunctor() : m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)),
    m_force(Registry::instance->configuration()->IBI_bins),
    m_potential(Registry::instance->configuration()->IBI_bins),
    m_exclusion_low(Registry::instance->configuration()->IBI_exclusion_low),
    m_exclusion_high(Registry::instance->configuration()->IBI_exclusion_high) { }

template<>
void IBI_ForceFunctor<IBI_Reload>::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("IBI force", pairList.size(),
        IBI_Force(pairList.getPairs(), pairList.getOffsets(), p_soa.f(), p_soa.r(),
            m_force.getXValues(), m_force.getYValues(), m_force.getLowerDefault(),
            m_force.getUpperDefault(), m_cutoff2, m_exclusion_low, m_exclusion_high));
}

template<>
void IBI_ForceFunctor<IBI_Reload>::IBI_Force::operator()(int idx) const {
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);
    const math::d3 shift0 = pair_offsets(idx, 0);
    const math::d3 shift1 = pair_offsets(idx, 1);

    const math::d3 r0 = r[s_idx_0] + shift0;
    const math::d3 r1 = r[s_idx_1] + shift1;

    if (math::pointInBox(r0, exclusion_low, exclusion_high) || math::pointInBox(r1, exclusion_low, exclusion_high)) return;

    const math::d3 dr = r0 - r1;
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;

    const auto r = Kokkos::sqrt(dr2);
    const auto F = FunctionPL::evaluateAt(r, force_def_low, force_def_high, force_x_values, force_y_values);
    const math::d3 force = (dr / r) * F;

    if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], force);
    if (shift1 == 0) Kokkos::atomic_sub(&f[s_idx_1], force);
}
