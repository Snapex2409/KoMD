//
// Created by alex on 7/31/24.
//

#include "Simulation.h"
#include "Registry.h"
#include "IO/Logging.h"

void Simulation::run() {
    auto config = Registry::instance->configuration();
    auto& integrators = Registry::instance->integrators();
    auto container = Registry::instance->moleculeContainer();
    auto& potentials = Registry::instance->forceFunctors();
    auto vtkWriter = Registry::instance->vtkWriter();

    const double dt = config->delta_t;
    const uint64_t write_freq = config->write_freq;
    vtkWriter->write("VTK_Output", 0);
    //initial pass to compute forces
    container->updateContainer();
    for (auto& potential : potentials) (*potential)();

    // main loop
    double t = 0;
    uint64_t simstep = 0;
    const uint64_t max_step = config->timesteps;
    while (simstep < max_step) {
        for (auto& integrator : integrators) integrator->integrate0();
        container->updateContainer();
        for (auto& potential : potentials) (*potential)();
        for (auto& integrator : integrators) integrator->integrate1();

        Log::simulation->info() << "Simstep=" << simstep << " t=" << t << std::endl;
        t += dt;
        simstep++;

        if (simstep % write_freq == 0) vtkWriter->write("VTK_Output", simstep);
    }
}
