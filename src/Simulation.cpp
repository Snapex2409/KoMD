//
// Created by alex on 7/31/24.
//

#include "Simulation.h"
#include "Registry.h"
#include "IO/Logging.h"
#include "IO/CheckpointIO.h"

void Simulation::run() {
    auto config = Registry::instance->configuration();
    auto& integrators = Registry::instance->integrators();
    auto& sensors = Registry::instance->sensors();
    auto temp_sens = Registry::instance->temperature_sensor();
    auto container = Registry::instance->moleculeContainer();
    auto& potentials = Registry::instance->forceFunctors();
    auto vtkWriter = Registry::instance->vtkWriter();
    auto& thermostats = Registry::instance->thermostats();

    const double dt = config->delta_t;
    const uint64_t write_freq = config->write_freq;
    vtkWriter->write("VTK_Output", 0);
    //initial pass to compute forces
    container->updateContainer();
    for (auto& potential : potentials) (*potential)();
    temp_sens->measure();
    Log::simulation->info() << "Initial temperature T=" << temp_sens->getTemperature() << std::endl;
    for (auto& thermostat : thermostats) thermostat->apply();

    // main loop
    double t = 0;
    uint64_t simstep = 0;
    const uint64_t max_step = config->timesteps;
    while (simstep < max_step) {
        for (auto& integrator : integrators) integrator->integrate0();
        container->updateContainer();
        for (auto& potential : potentials) (*potential)();
        for (auto& integrator : integrators) integrator->integrate1();

        for (auto& sensor : sensors) sensor->measure();

        Log::simulation->info() << "Simstep=" << simstep << " t=" << t << " T=" << temp_sens->getTemperature() << std::endl;
        t += dt;
        simstep++;

        if (simstep % write_freq == 0) vtkWriter->write("VTK_Output", simstep);

        for (auto& thermostat : thermostats) thermostat->apply();
    }

    for (auto& sensor : sensors) sensor->write(simstep);
    if (config->storeCheckpoint) CheckpointIO::writeCheckpoint(simstep);
}
