//
// Created by alex on 8/4/24.
//

#include "LJ12_6.h"
#include "Registry.h"

LJ12_6::LJ12_6() : m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)) { }

void LJ12_6::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("LJ12-6", pairList.size(),
                         LJ12_6_Force(pairList.getPairs(), pairList.getOffsets(), p_soa.f(), p_soa.id(), p_soa.r(), p_soa.sigma(), p_soa.epsilon(), m_cutoff2));
}

void LJ12_6::LJ12_6_Force::operator()(int idx) const {
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);
    
    if (id(s_idx_0) == id(s_idx_1)) return;
    const math::d3 shift0 = pair_offsets(idx, 0);
    const math::d3 shift1 = pair_offsets(idx, 1);

    const math::d3 dr = (r[s_idx_0] + shift0) - (r[s_idx_1] + shift1);
    const double dr2 = dr.dot(dr);
    const int id0 = id[s_idx_0];
    const int id1 = id[s_idx_1];
    //printf("%d %d cut %d, dr: %f, shift 0:%f %f %f, shift 1:%f %f %f\n", id0, id1, dr2 > cutoff2, dr2, shift0.x(), shift0.y(), shift0.z(), shift1.x(), shift1.y(), shift1.z());
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps[s_idx_0] * eps[s_idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const auto force = dr * fac;

    if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], force);
    if (shift1 == 0) Kokkos::atomic_sub(&f[s_idx_1], force);
    //printf("%d %d fx %f fy %f\n", id0, id1, dr.x()*fac, dr.y()*fac);
}
