//
// Created by alex on 8/21/24.
//

#include "OCLJ.h"
#include "Registry.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

OCLJ::OCLJ() : m_cutoff2(std::pow(Registry::instance->configuration()->cutoff, 2)) {
    p_run_pairs = false;
    auto config = Registry::instance->configuration();
    m_domain_size = config->domainHigh - config->domainLow;
    m_container = std::dynamic_pointer_cast<OneCell>(Registry::instance->moleculeContainer());
    m_low_bound = config->domainLow + config->cutoff;
}

void OCLJ::handleCell(Cell &cell) {
    Kokkos::parallel_for("OCLJ - Cell", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
            {0,0,0},
            {static_cast<long>(cell.getNumIndices()), static_cast<long>(cell.getNumIndices()), 8}),
            OCLJ_Force(p_soa.f(), p_soa.id(), p_soa.r(), m_container->getCOM(), p_soa.sigma(), p_soa.epsilon(),
                       m_cutoff2, m_domain_size, m_low_bound));
}

void OCLJ::OCLJ_Force::operator()(int idx_0, int idx_1, int s_idx) const {
    if (idx_0 == idx_1) return; // do not compute for same site, we do not use Newton 3 here
    if (id(idx_0) == id(idx_1)) return; // must be different molecule

    const math::i3 shift_dir {s_idx & 0b1, (s_idx & 0b10) >> 1, (s_idx & 0b100) >> 2};
    const math::d3 shift = (com[idx_0] <= low_bound) * shift_dir * domain_size;
    if (s_idx > 0 && shift == 0) return; // no shift performed, when shift was intended

    const math::d3 dr = (r[idx_0] + shift) - r[idx_1];
    const double dr2 = dr.dot(dr);
    if (dr2 > cutoff2) return;
    const double invdr2 = 1. / dr2;

    const double sigma = (sig[idx_0] + sig[idx_1]) / 2.0;
    const double epsilon = Kokkos::sqrt(eps[idx_0] * eps[idx_1]);

    const double sig2 = sigma * sigma;
    const double lj6 = Kokkos::pow(sig2 * invdr2, 3.0);
    const double lj12 = Kokkos::pow(lj6, 2.0);
    const double fac = 24.0 * epsilon * (2.0 * lj12 - lj6) * invdr2;
    const math::d3 force = dr * fac;
    Kokkos::atomic_add(&f[idx_0], force);
    if (s_idx > 0) {
        Kokkos::atomic_sub(&f[idx_1], force);
    }
}
