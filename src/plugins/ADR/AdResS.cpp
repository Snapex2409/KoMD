//
// Created by alex on 10/15/24.
//

#include "AdResS.h"

#include "Registry.h"
#include "math/Geometry.h"
#include "ADR_Force.h"
#include "ADR_Integrator.h"

AdResS::AdResS() : Plugin("ADR"), m_fp_low(Registry::instance->configuration()->ADR_low),
                   m_fp_high(Registry::instance->configuration()->ADR_high),
                   m_h_dim(Registry::instance->configuration()->ADR_h_dim), m_total_mass(0), m_site_count(1), m_is_first(true) { }

void AdResS::init() {
    // check components for validity
    auto& components = Registry::instance->components();
    if (components.size() > 1) throw std::runtime_error("For ADR only 1 Hybrid component allowed.");
    if (components[0].getSites().size() < 2) throw std::runtime_error("We require at least one CG and one FP site");
    if (components[0].getSites()[0].getMass() != 0) throw std::runtime_error("First site must be mass-less CG interaction site.");
    for (int idx = 1; idx < components[0].getSites().size(); idx++) {
        auto& site = components[0].getSites()[idx];
        if (site.getMass() == 0) throw std::runtime_error("Only one mass-less site allowed");
        m_total_mass += site.getMass();
        m_site_count++;
    }
    m_weights = KW::vec_t<double>("ADR - weights", Registry::instance->moleculeContainer()->getSOA().size());

    Registry::instance->forceFunctors().clear();
    Registry::instance->forceFunctors().push_back(std::make_unique<ADR_Force>(m_weights));

    Registry::instance->integrators().clear();
    Registry::instance->integrators().push_back(std::make_unique<ADR_Integrator>());
}

void AdResS::pre_main_loop() {
    m_is_first = false;
}

void AdResS::pre_container_update() {
    if (m_is_first) return;

    auto container = Registry::instance->moleculeContainer();
    const auto coms = container->getCOM();
    auto& soa = container->getSOA();
    Kokkos::parallel_for("ADR - CG Pos", soa.size(), CG_Pos_Kernel(soa.r(), coms, soa.mass()));
    Kokkos::fence("ADR - CG Pos fence");
}

void AdResS::post_container_update() {
    const auto coms = Registry::instance->moleculeContainer()->getCOM();
    Kokkos::parallel_for("ADR - weight calc", coms.size(), Weight_Kernel(coms, m_weights, m_fp_low, m_fp_high, m_h_dim));
    Kokkos::fence("ADR - weight fence");
}

void AdResS::post_forces() {
    auto& soa = Registry::instance->moleculeContainer()->getSOA();
    Kokkos::parallel_for("ADR - force distribution", soa.size(), Force_Distribute_Kernel(soa.f(), soa.mass(), m_site_count, m_total_mass));
    Kokkos::fence("ADR - force distribution fence");
}

double AdResS::weight(const math::d3 &r, const math::d3 &fp_low, const math::d3 &fp_high, const math::d3 &h_dim) {
    if (math::pointInBox(r, fp_low, fp_high)) return 1.0;   // is FP
    if (!math::pointInBox(r, fp_low - h_dim, fp_high + h_dim)) return 0.0; // is CG

    // Handle H
    const math::d3 dist_dir = math::max(
        math::max(fp_low - r, {0, 0, 0}),
        math::max(r - fp_high, {0, 0, 0})
    );
    const double dist = dist_dir.L2();
    const math::d3 active_h_mask = math::d3{1, 1, 1} - (dist_dir <= math::d3{0, 0, 0}) * (dist_dir >= math::d3{0, 0, 0});
    const double hDim = math::pow(active_h_mask * h_dim, 2.0).L2();
    if (dist >= hDim) return 0.0; // outside of rounded region -> treat as CG
    return Kokkos::pow(Kokkos::cos(M_PI/(2*hDim) * dist), 2.0);
}

void AdResS::Weight_Kernel::operator()(int idx) const {
    weights(idx) = AdResS::weight(coms(idx), fp_low, fp_high, h_dim);
}

void AdResS::Force_Distribute_Kernel::operator()(int idx) const {
    const double mass = m(idx);
    if (mass == 0) return;

    const int offset = idx % site_count;
    const double frac = mass / total_mass;

    const auto force = f(idx - offset) * frac;
    f(idx) += force;
}

void AdResS::CG_Pos_Kernel::operator()(int idx) const {
    if (m(idx) != 0) return;
    r(idx) = com(idx);
}
