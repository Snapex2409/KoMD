//
// Created by alex on 8/4/24.
//

#include "ATM2B.h"
#include "Registry.h"
#include "util/constants.h"

ATM2B::ATM2B() : m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)) {
    auto config = Registry::instance->configuration();
    double nu = config->energy_3b;
    nu = nu * 1e+9 * Constants::conv_J_Ei; // convert external unit to internal
    const double rho = Registry::instance->moleculeContainer()->getNumMolecules() / (config->domainHigh - config->domainLow).product();
    // we assume we only have 1 component system
    const double eps = Registry::instance->components()[0].getSites()[0].getEpsilon() * Constants::conv_J_Ei;
    const double sig = Registry::instance->components()[0].getSites()[0].getSigma();
    const double sig6 = std::pow(sig, 6);
    m_potential_factor = 1.0 - ((2*nu*rho)/(3*eps*sig6));
}

void ATM2B::handlePairList(PairList &pairList) {
    Kokkos::parallel_for("ATM2B", pairList.size(),
                         ATM2B_Force(pairList.getPairs(), pairList.getOffsets(), p_soa.f(), p_soa.id(), p_soa.r(),
                                   p_soa.sigma(), p_soa.epsilon(), m_cutoff2, m_potential_factor));
}

void ATM2B::ATM2B_Force::operator()(int idx) const {
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
    const auto force = dr * fac * potential_factor;

    if (shift0 == 0) Kokkos::atomic_add(&f[s_idx_0], force);
    if (shift1 == 0) Kokkos::atomic_sub(&f[s_idx_1], force);
    //printf("%d %d fx %f fy %f\n", id0, id1, dr.x()*fac, dr.y()*fac);
}
