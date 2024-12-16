//
// Created by alex on 7/31/24.
//

#ifndef KOMD_SIMULATION_H
#define KOMD_SIMULATION_H

#include <memory>

#include "IO/Configuration.h"

/**
 * Contains main simulation loop, that calls all other components.
 * */
class Simulation {
public:
    /**
     * Runs the simulation
     * */
    void run();

    /**
     * Stop the simulation at the end of this timestep.
     * */
    void stop();

    uint64_t simstep() const { return m_simstep; }

private:
    /// maximum time step
    uint64_t m_max_step = 0;
    /// current time step
    uint64_t m_simstep = 0;
};


#endif //KOMD_SIMULATION_H
