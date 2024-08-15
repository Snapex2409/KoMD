//
// Created by alex on 7/31/24.
//

#include "Integrator.h"
#include "Registry.h"
#include "molecule/Molecule.h"

Integrator::Integrator() : m_delta_t(Registry::instance->configuration()->delta_t) { }

void Integrator::integrate0() {
    auto container = Registry::instance->moleculeContainer();

    const double dt_halve = m_delta_t * 0.5;
    for (auto it = container->iterator(MoleculeContainer::SITE, MoleculeContainer::DOMAIN); it.isValid(); ++it) {
        const double dtInv2m = dt_halve / it.mass();
        it.v() += it.f() * dtInv2m;
        it.r() += it.v() * m_delta_t;
    }
}

void Integrator::integrate1() {
    auto container = Registry::instance->moleculeContainer();

    const double dt_halve = m_delta_t * 0.5;
    for (auto it = container->iterator(MoleculeContainer::SITE, MoleculeContainer::DOMAIN); it.isValid(); ++it) {
        const double dtInv2m = dt_halve / it.mass();
        it.v() += it.f() * dtInv2m;
    }
}
