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

std::vector<std::shared_ptr<Sensor>> &Registry::sensors() {
    return m_sensors;
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

std::shared_ptr<VTKWriter> Registry::vtkWriter() {
    return m_vtk_writer;
}

std::shared_ptr<VTKWriter> &Registry::vtkWriter_ptr() {
    return m_vtk_writer;
}

std::shared_ptr<LJ12_6_Sensor> Registry::potential_sensor() {
    return m_sensor_pot;
}

std::shared_ptr<LJ12_6_Sensor> & Registry::potential_sensor_ptr() {
    return m_sensor_pot;
}

std::shared_ptr<TemperatureSensor> Registry::temperature_sensor() {
    return m_sensor_temp;
}

std::shared_ptr<TemperatureSensor> & Registry::temperature_sensor_ptr() {
    return m_sensor_temp;
}

std::vector<std::unique_ptr<Thermostat>> & Registry::thermostats() {
    return m_thermostats;
}

std::shared_ptr<DisplacementSensor> Registry::displacement_sensor() {
    return m_sensor_disp;
}

std::shared_ptr<DisplacementSensor> &Registry::displacement_sensor_ptr() {
    return m_sensor_disp;
}
