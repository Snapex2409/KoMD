//
// Created by alex on 10/15/24.
//

#include "ADR_Integrator.h"
#include "Registry.h"

ADR_Integrator::ADR_Integrator(double total_mass) : m_total_mass(total_mass) { }

void ADR_Integrator::integrate0() {
    auto container = Registry::instance->moleculeContainer();
    SOA& soa = container->getSOA();

    const double dt_halve = p_delta_t * 0.5;
    Kokkos::parallel_for("ADR Integrate 0", soa.size(), ADR_Step0(soa.r(), soa.v(), soa.f(), soa.mass(), dt_halve, p_delta_t, m_total_mass));
    Kokkos::fence("ADR Integration0 fence");
}

void ADR_Integrator::integrate1() {
    auto container = Registry::instance->moleculeContainer();

    const double dt_halve = p_delta_t * 0.5;
    SOA& soa = container->getSOA();
    Kokkos::parallel_for("ADR Integrate 1", soa.size(), ADR_Step1(soa.v(), soa.f(), soa.mass(), dt_halve, m_total_mass));
    Kokkos::fence("ADR Integration1 fence");
}

void ADR_Integrator::ADR_Step0::operator()(int idx) const {
    const double mass = m(idx);
    double active_mass = mass;
    if (mass == 0) active_mass = total_mass;

    const double dtInv2m = dt_halve / active_mass;
    v(idx) += f(idx) * dtInv2m;

    if (mass == 0) return; // we do not integrate CG-r here
    r(idx) += v(idx) * dt;
}

void ADR_Integrator::ADR_Step1::operator()(int idx) const {
    const double mass = m(idx);
    double active_mass = mass;
    if (mass == 0) active_mass = total_mass;

    const double dtInv2m = dt_halve / active_mass;
    v(idx) += f(idx) * dtInv2m;
}
