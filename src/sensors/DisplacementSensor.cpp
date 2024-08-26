//
// Created by alex on 8/22/24.
//

#include "DisplacementSensor.h"
#include "Registry.h"

DisplacementSensor::DisplacementSensor() : Sensor("Displacement"), m_lambda(1) {
    auto config = Registry::instance->configuration();
    const double density = config->density;
    const math::d3 domain_size = config->domainHigh - config->domainLow;
    const double domain_volume = domain_size.product();
    const double num_molecules = density * domain_volume;
    const double mol_per_dim = std::pow(num_molecules, 1./3.);
    const math::d3 spacing = domain_size / mol_per_dim;
    m_lattice_unit_size = spacing;

    m_com = KW::vec_t<math::d3>("Displacement Sensor CoM", Registry::instance->moleculeContainer()->getNumMolecules());
}

void DisplacementSensor::measure() {
    auto container = Registry::instance->moleculeContainer();
    const auto count = container->getNumMolecules();
    container->getCenterOfMassPositions(m_com);

    double lx = 0, ly = 0, lz = 0;
    Kokkos::parallel_reduce("Displacement Sensor", count, Displacement_Kernel(m_com, m_lattice_unit_size), lx, ly, lz);
    Kokkos::fence("Displacement fence");

    m_lambda = 1.0 / (3.0 * count) * (lx + ly + lz);
}

void DisplacementSensor::write(uint64_t simstep) { }

void DisplacementSensor::Displacement_Kernel::operator()(int idx, double& lx, double& ly, double& lz) const {
    const math::d3 r = com[idx];
    lx += 0.5 * Kokkos::cos(2.0 * M_PI * (r.x() - 0.5) / lattice_unit_size.x()) + 0.5;
    ly += 0.5 * Kokkos::cos(2.0 * M_PI * (r.y() - 0.5) / lattice_unit_size.y()) + 0.5;
    lz += 0.5 * Kokkos::cos(2.0 * M_PI * (r.z() - 0.5) / lattice_unit_size.z()) + 0.5;
}
