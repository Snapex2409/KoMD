//
// Created by alex on 7/31/24.
//

#include "Integrator.h"
#include "Registry.h"

Integrator::Integrator() : p_delta_t(Registry::instance->configuration()->delta_t) { }

void Integrator::integrate0() {
    auto container = Registry::instance->moleculeContainer();
    SOA& soa = container->getSOA();

    const double dt_halve = p_delta_t * 0.5;
    Kokkos::parallel_for("Integrate 0", soa.size(), Step0(soa.r(), soa.v(), soa.f(), soa.mass(), dt_halve, p_delta_t));
    Kokkos::fence("Integration0 fence");
}

void Integrator::integrate1() {
    auto container = Registry::instance->moleculeContainer();

    const double dt_halve = p_delta_t * 0.5;
    SOA& soa = container->getSOA();
    Kokkos::parallel_for("Integrate 1", soa.size(), Step1(soa.v(), soa.f(), soa.mass(), dt_halve));
    Kokkos::fence("Integration1 fence");
}

void Integrator::Step0::operator()(int idx) const {
    const double dtInv2m = dt_halve / m(idx);
    v(idx) += f(idx) * dtInv2m;
    r(idx) += v(idx) * dt;
}

void Integrator::Step1::operator()(int idx) const {
    const double dtInv2m = dt_halve / m(idx);
    v(idx) += f(idx) * dtInv2m;
}
