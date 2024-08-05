//
// Created by alex on 7/31/24.
//

#include "Registry.h"


std::unique_ptr<Registry> Registry::instance;

Registry::Registry() : m_configuration(), m_forceFunctors(), m_integrators(), m_moleculeContainer(), m_simulation() {
    if (instance != nullptr) throw std::runtime_error("Registry already exists");
}

std::shared_ptr<Configuration> Registry::configuration() {
    return m_configuration;
}

std::vector<std::unique_ptr<ForceFunctor>> &Registry::forceFunctors() {
    return m_forceFunctors;
}

std::vector<std::unique_ptr<Integrator>> &Registry::integrators() {
    return m_integrators;
}

std::shared_ptr<MoleculeContainer> Registry::moleculeContainer() {
    return m_moleculeContainer;
}

std::shared_ptr<Simulation> Registry::simulation() {
    return m_simulation;
}

std::shared_ptr<Configuration> &Registry::configuration_ptr() {
    return m_configuration;
}

std::shared_ptr<MoleculeContainer> &Registry::moleculeContainer_ptr() {
    return m_moleculeContainer;
}

std::shared_ptr<Simulation> &Registry::simulation_ptr() {
    return m_simulation;
}

std::shared_ptr<Boundary> Registry::boundary() {
    return m_boundary;
}

std::shared_ptr<Boundary> &Registry::boundary_ptr() {
    return m_boundary;
}