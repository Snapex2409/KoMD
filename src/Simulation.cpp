//
// Created by alex on 7/31/24.
//

#include "Simulation.h"
#include "Registry.h"

void Simulation::run() {
    auto config = Registry::instance->configuration();
    const double dt = config->delta_t;
    //initial pass to compute forces

    // main loop
    double t = 0;
    uint64_t simstep = 0;
    const uint64_t max_step = config->timesteps;
    while (simstep < max_step) {


        t += dt;
        simstep++;
    }
}
