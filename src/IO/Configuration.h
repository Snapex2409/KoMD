//
// Created by alex on 7/31/24.
//

#ifndef KOMD_CONFIGURATION_H
#define KOMD_CONFIGURATION_H

#include "util/defaults.h"
#include "math/Array.h"

struct Configuration {
    double temperature = Defaults::temperature;
    double delta_t = Defaults::delta_t;
    double cutoff = Defaults::cutoff;
    double density = Defaults::density;
    double epsilon = Defaults::epsilon;
    double sigma = Defaults::sigma;
    double stiffness_factor = Defaults::stiffness_factor;
    double limit_factor = Defaults::limit_factor;
    math::d3 domainLow = Defaults::domainLow;
    math::d3 domainHigh = Defaults::domainHigh;
    bool enableSOA = Defaults::enableSOA;
    uint64_t timesteps = Defaults::timesteps;
    uint64_t write_freq = Defaults::write_freq;
};


#endif //KOMD_CONFIGURATION_H
