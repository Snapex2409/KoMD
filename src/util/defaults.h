//
// Created by alex on 7/31/24.
//

#ifndef KOMD_DEFAULTS_H
#define KOMD_DEFAULTS_H

#include "math/Array.h"

namespace Defaults {
    static constexpr double sigma = 1.0;                // 1 Armstrong
    static constexpr double epsilon = 1.0;              // 1 Joule
    static constexpr double mass = 1.0;                 // 1 Dalton = 1u
    static constexpr double delta_t = 0.001;            // 1 ps, 0.001 = 1 fs
    static constexpr double stiffness_factor = 30;      // unit-less
    static constexpr double density = 0.1;              // unit-less, number density
    static constexpr math::d3 r {0, 0, 0};
    static constexpr math::d3 v {0, 0, 0};
    static constexpr double temperature = 300;          // 300 Kelvin
    static constexpr math::d3 domainLow {0, 0, 0};      // in Armstrong
    static constexpr math::d3 domainHigh {10, 10, 10};  // in Armstrong
    static constexpr double cutoff = sigma * 2.5;       // in Armstrong
    static constexpr bool enableSOA = false;
}

#endif //KOMD_DEFAULTS_H
