//
// Created by alex on 7/31/24.
//

#include "Integrator.h"
#include "Registry.h"

#include "Kokkos_Core.hpp"

Integrator::Integrator() : m_delta_t(Registry::instance->configuration()->delta_t) { }

void Integrator::integrate0() {
    auto container = Registry::instance->moleculeContainer();

    const double dt_halve = m_delta_t * 0.5;
    for (auto it = container->iteratorCell(MoleculeContainer::DOMAIN); it.isValid(); ++it) {
        Cell& cell = it.cell();
        SOA& soa = cell.soa();
        Kokkos::parallel_for("Integrate 0", soa.size(), Step0(soa, dt_halve, m_delta_t));
    }
    Kokkos::fence("Integration0 fence");
}

void Integrator::integrate1() {
    auto container = Registry::instance->moleculeContainer();

    const double dt_halve = m_delta_t * 0.5;
    for (auto it = container->iteratorCell(MoleculeContainer::DOMAIN); it.isValid(); ++it) {
        Cell& cell = it.cell();
        SOA& soa = cell.soa();
        Kokkos::parallel_for("Integrate 1", soa.size(), Step1(soa, dt_halve));
    }
    Kokkos::fence("Integration1 fence");
}

void Integrator::Step0::operator()(int idx) const {
    auto& v = soa.v();
    auto& f = soa.f();
    auto& r = soa.r();
    auto& m = soa.mass();

    const double dtInv2m = dt_halve / m(idx);
    v(idx) += f(idx) * dtInv2m;
    r(idx) += v(idx) * dt;
}

void Integrator::Step1::operator()(int idx) const {
    auto& v = soa.v();
    auto& f = soa.f();
    auto& m = soa.mass();

    const double dtInv2m = dt_halve / m(idx);
    v(idx) += f(idx) * dtInv2m;
}
