//
// Created by alex on 7/31/24.
//

#ifndef KOMD_REGISTRY_H
#define KOMD_REGISTRY_H

#include <memory>
#include <vector>
#include <thermostats/Thermostat.h>

#include "container/MoleculeContainer.h"
#include "integrators/Integrator.h"
#include "IO/Configuration.h"
#include "potentials/ForceFunctor.h"
#include "Simulation.h"
#include "boundary/Boundary.h"
#include "IO/VTKWriter.h"
#include "sensors/Sensor.h"
#include "sensors/TemperatureSensor.h"

class Registry {
public:
    Registry();
    static std::unique_ptr<Registry> instance;

    std::shared_ptr<Configuration> configuration();
    std::vector<std::unique_ptr<ForceFunctor>>& forceFunctors();
    std::vector<std::unique_ptr<Integrator>>& integrators();
    std::vector<std::shared_ptr<Sensor>>& sensors();
    std::shared_ptr<MoleculeContainer> moleculeContainer();
    std::shared_ptr<Simulation> simulation();
    std::shared_ptr<Boundary> boundary();
    std::shared_ptr<VTKWriter> vtkWriter();
    std::shared_ptr<TemperatureSensor> temperature_sensor();
    std::vector<std::unique_ptr<Thermostat>>& thermostats();

    std::shared_ptr<Configuration>& configuration_ptr();
    std::shared_ptr<MoleculeContainer>& moleculeContainer_ptr();
    std::shared_ptr<Simulation>& simulation_ptr();
    std::shared_ptr<Boundary>& boundary_ptr();
    std::shared_ptr<VTKWriter>& vtkWriter_ptr();
    std::shared_ptr<TemperatureSensor>& temperature_sensor_ptr();
private:
    std::shared_ptr<Configuration> m_configuration;
    std::vector<std::unique_ptr<ForceFunctor>> m_forceFunctors;
    std::vector<std::unique_ptr<Integrator>> m_integrators;
    std::vector<std::shared_ptr<Sensor>> m_sensors;
    std::shared_ptr<MoleculeContainer> m_moleculeContainer;
    std::shared_ptr<Simulation> m_simulation;
    std::shared_ptr<Boundary> m_boundary;
    std::shared_ptr<VTKWriter> m_vtk_writer;
    std::shared_ptr<TemperatureSensor> m_sensor_temp;
    std::vector<std::unique_ptr<Thermostat>> m_thermostats;
};


#endif //KOMD_REGISTRY_H
