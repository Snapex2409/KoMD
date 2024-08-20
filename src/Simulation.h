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
};


#endif //KOMD_SIMULATION_H
