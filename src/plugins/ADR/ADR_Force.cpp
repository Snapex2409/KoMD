//
// Created by alex on 10/15/24.
//

#include "ADR_Force.h"
#include "Registry.h"

ADR_Force::ADR_Force(KW::vec_t<double>& weights) : ForceFunctor(),
    m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)),
    m_stiffness_factor(Registry::instance->configuration()->stiffness_factor),
    m_weights_ref(weights) { }

void ADR_Force::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("ADR - Force", pairList.size(),
        Force_Kernel(pairList.getPairs(), pairList.getOffsets(), p_soa.f(), p_soa.id(), p_soa.r(),
            p_soa.sigma(), p_soa.epsilon(), p_soa.mass(), m_cutoff2, m_stiffness_factor, m_weights_ref));
}

void ADR_Force::Force_Kernel::operator()(int idx) const {
    const uint64_t s_idx_0 = pairs(idx, 0);
    const uint64_t s_idx_1 = pairs(idx, 1);

    if (id(s_idx_0) == id(s_idx_1)) return handle_FENE(idx, s_idx_0, s_idx_1);
    return handle_LJ(idx, s_idx_0, s_idx_1);
}

void ADR_Force::Force_Kernel::handle_FENE(int idx, uint64_t s_idx_0, uint64_t s_idx_1) const {
    const int res0 = mass(s_idx_0) != 0; // res 1 is FP
    const int res1 = mass(s_idx_1) != 0;
    if (res0 != res1) return;
    if (res0 == 0) return;
    if (pair_offsets(idx, 0) != 0 || pair_offsets(idx, 1) != 0) return;
    const double w = weights(s_idx_0);
    if (w == 0) return; // no FENE for CG

    const math::d3 dr = r[s_idx_1] - r[s_idx_0];
    const math::d3 invdr = math::d3{1, 1, 1} / dr;
    const double dr2 = dr.dot(dr);

    const double sigma = (sig[s_idx_0] + sig[s_idx_1]) / 2.0;
    const double R02 = std::pow(1.5 * sigma, 2);
    if (R02 <= dr2) {
        const auto force = invdr * 1e+30;
        Kokkos::atomic_add(&f[s_idx_0], force);
        Kokkos::atomic_sub(&f[s_idx_1], force);
        return;
    }

    const double epsilon = std::sqrt(eps[s_idx_0] * eps[s_idx_1]);
    const double k = stiffness_factor * epsilon / std::pow(sigma, 2);

    const double fac = k / (1 - (dr2/R02));
    const auto force = dr * fac;

    Kokkos::atomic_add(&f[s_idx_0], force);
    Kokkos::atomic_sub(&f[s_idx_1], force);
}

void ADR_Force::Force_Kernel::handle_LJ(int idx, uint64_t s_idx_0, uint64_t s_idx_1) const {
    const int res0 = mass(s_idx_0) != 0;
    const int res1 = mass(s_idx_1) != 0;
    if (res0 != res1) return;

    const math::d3 shift0 = pair_offsets(idx, 0);
    const math::d3 shift1 = pair_offsets(idx, 1);

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

    const double w0 = weights(s_idx_0);
    const double w1 = weights(s_idx_1);
    const double w = res0 * w0 * w1 + (1 - res0) * (1 - w0 * w1);

    const auto force = dr * (fac * w);

    if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], force);
    if (shift1 == 0) Kokkos::atomic_sub(&f[s_idx_1], force);
}
