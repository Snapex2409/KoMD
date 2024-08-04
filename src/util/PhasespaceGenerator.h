//
// Created by alex on 8/4/24.
//

#ifndef KOMD_PHASESPACEGENERATOR_H
#define KOMD_PHASESPACEGENERATOR_H


#include "math/Array.h"

class PhasespaceGenerator {
public:
    static void generate();
private:
    /**
     * Generates a random sample of the Maxwell Boltzmann distribution.
     * @param T temperature in Kelvin
     * @param m mass in Dalton (u)
     * */
    static math::d3 getMaxwellBoltzmannVelocity(double T, double m);
};


#endif //KOMD_PHASESPACEGENERATOR_H
