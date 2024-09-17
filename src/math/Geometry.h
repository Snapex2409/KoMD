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
    KOKKOS_INLINE_FUNCTION bool pointInBox(const d3& point, const d3& low, const d3& high) {
        return point.x() >= low.x() && point.y() >= low.y() && point.z() >= low.z() &&
               point.x() <= high.x() && point.y() <= high.y() && point.z() <= high.z();
    }

    /**
     * Checks if the boxes defines by the pairs (low_0, high_0) and (low_1, high_1) have any sort of intersection
     * */
    KOKKOS_INLINE_FUNCTION bool boxesIntersect(const d3& low_0, const d3& high_0, const d3& low_1, const d3& high_1) {
        return high_0.x() >= low_1.x() && high_1.x() >= low_0.x() &&
               high_0.y() >= low_1.y() && high_1.y() >= low_0.y() &&
               high_0.z() >= low_1.z() && high_1.z() >= low_0.z();
    }

    /**
     * Cantor mapping function from N² to N
     * */
    template<typename T>
    KOKKOS_INLINE_FUNCTION T cantor2(T x, T y) {
        return (x + y) * (x + y + 1) / 2 + y;
    }

    /**
     * Cantor mapping function for N³ to N
     * */
    template<typename T>
     KOKKOS_INLINE_FUNCTION T cantor3(const Array<T, 3>& coord) {
         return cantor2(cantor2(coord.x(), coord.y()), coord.z());
     }
}

#endif //KOMD_GEOMETRY_H
