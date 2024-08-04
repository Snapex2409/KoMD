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
}
#endif //KOMD_CONSTANTS_H
