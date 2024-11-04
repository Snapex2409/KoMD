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
#include "potentials/ForceFunctor3B.h"
#include "Simulation.h"
#include "IO/VTKWriter.h"
#include "sensors/Sensor.h"
#include "sensors/TemperatureSensor.h"
#include "sensors/LJ12_6_Sensor.h"
#include "thermostats/Thermostat.h"
#include "sensors/DisplacementSensor.h"
#include "molecule/Component.h"
#include "plugins/Plugin.h"

class Registry {
public:
    Registry();
    static std::unique_ptr<Registry> instance;

    std::shared_ptr<Configuration> configuration();
    std::vector<std::unique_ptr<ForceFunctor>>& forceFunctors();
    std::vector<std::unique_ptr<ForceFunctor3B>>& forceFunctors3b();
    std::vector<std::unique_ptr<Integrator>>& integrators();
    std::vector<std::shared_ptr<Sensor>>& sensors();
    std::shared_ptr<MoleculeContainer> moleculeContainer();
    std::shared_ptr<Simulation> simulation();
    std::shared_ptr<VTKWriter> vtkWriter();
    std::shared_ptr<TemperatureSensor> temperature_sensor();
    std::shared_ptr<LJ12_6_Sensor> potential_sensor();
    std::shared_ptr<DisplacementSensor> displacement_sensor();
    std::vector<std::unique_ptr<Thermostat>>& thermostats();
    std::vector<Component>& components();
    std::vector<std::unique_ptr<Plugin>>& plugins();

    std::shared_ptr<Configuration>& configuration_ptr();
    std::shared_ptr<MoleculeContainer>& moleculeContainer_ptr();
    std::shared_ptr<Simulation>& simulation_ptr();
    std::shared_ptr<VTKWriter>& vtkWriter_ptr();
    std::shared_ptr<TemperatureSensor>& temperature_sensor_ptr();
    std::shared_ptr<LJ12_6_Sensor>& potential_sensor_ptr();
    std::shared_ptr<DisplacementSensor>& displacement_sensor_ptr();
private:
    std::shared_ptr<Configuration> m_configuration;
    std::vector<std::unique_ptr<ForceFunctor>> m_forceFunctors;
    std::vector<std::unique_ptr<ForceFunctor3B>> m_forceFunctors3b;
    std::vector<std::unique_ptr<Integrator>> m_integrators;
    std::vector<std::shared_ptr<Sensor>> m_sensors;
    std::shared_ptr<MoleculeContainer> m_moleculeContainer;
    std::shared_ptr<Simulation> m_simulation;
    std::shared_ptr<VTKWriter> m_vtk_writer;
    std::shared_ptr<TemperatureSensor> m_sensor_temp;
    std::shared_ptr<LJ12_6_Sensor> m_sensor_pot;
    std::shared_ptr<DisplacementSensor> m_sensor_disp;
    std::vector<std::unique_ptr<Thermostat>> m_thermostats;
    std::vector<Component> m_components;
    std::vector<std::unique_ptr<Plugin>> m_plugins;
};


#endif //KOMD_REGISTRY_H
