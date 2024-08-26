//
// Created by alex on 8/4/24.
//

#ifndef KOMD_CONSTANTS_H
#define KOMD_CONSTANTS_H

#include "math/SciValue.h"

namespace Constants {
    /// boltzmann constant
    static constexpr SciValue kB = SciValue(1.380649, -23);
    /// dalton to kg factor f: 1 Da = f * kg
    static constexpr SciValue conv_Da_kg = SciValue(1.66053906892, -27);
    /// m/s to A/ps
    static constexpr SciValue conv_ms_Aps = SciValue(1.0, -2);
    /// A/ps to m/s
    static constexpr SciValue conv_Aps_ms = SciValue(1.0, 2);
    /// Joule to internal energy (u*A²/ps²)
    static constexpr SciValue conv_J_Ei = SciValue(6.022140754, 22);
}
#endif //KOMD_CONSTANTS_H
