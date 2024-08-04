//
// Created by alex on 8/2/24.
//

#ifndef KOMD_GEOMETRY_H
#define KOMD_GEOMETRY_H

#include "Array.h"

namespace math {
    /**
     * Checks if the specified point is within the box defined by low and high.
     * */
    bool pointInBox(const d3& point, const d3& low, const d3& high);
}

#endif //KOMD_GEOMETRY_H
