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

    bool boxesIntersect(const d3& low_0, const d3& high_0, const d3& low_1, const d3& high_1);
}

#endif //KOMD_GEOMETRY_H
