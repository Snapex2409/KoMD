//
// Created by alex on 7/31/24.
//

#ifndef KOMD_REGISTRY_H
#define KOMD_REGISTRY_H

#include <memory>
#include <vector>

#include "container/MoleculeContainer.h"
#include "integrators/Integrator.h"
#include "IO/Configuration.h"
#include "potentials/ForceFunctor.h"
#include "Simulation.h"
#include "boundary/Boundary.h"

class Registry {
public:
    Registry();
    static std::unique_ptr<Registry> instance;

    std::shared_ptr<Configuration> configuration();
    std::vector<ForceFunctor>& forceFunctors();
    std::vector<Integrator>& integrators();
    std::shared_ptr<MoleculeContainer> moleculeContainer();
    std::shared_ptr<Simulation> simulation();
    std::shared_ptr<Boundary> boundary();

    std::shared_ptr<Configuration>& configuration_ptr();
    std::shared_ptr<MoleculeContainer>& moleculeContainer_ptr();
    std::shared_ptr<Simulation>& simulation_ptr();
    std::shared_ptr<Boundary>& boundary_ptr();
private:
    std::shared_ptr<Configuration> m_configuration;
    std::vector<ForceFunctor> m_forceFunctors;
    std::vector<Integrator> m_integrators;
    std::shared_ptr<MoleculeContainer> m_moleculeContainer;
    std::shared_ptr<Simulation> m_simulation;
    std::shared_ptr<Boundary> m_boundary;
};


#endif //KOMD_REGISTRY_H
