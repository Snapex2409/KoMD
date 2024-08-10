//
// Created by alex on 8/2/24.
//
#include "Geometry.h"

bool math::pointInBox(const math::d3 &point, const math::d3 &low, const math::d3 &high) {
    return point.x() >= low.x() && point.y() >= low.y() && point.z() >= low.z() &&
           point.x() <= high.x() && point.y() <= high.y() && point.z() <= high.z();
}

bool math::boxesIntersect(const d3 &low_0, const d3 &high_0, const d3 &low_1, const d3 &high_1) {
    return high_0.x() >= low_1.x() && high_1.x() >= low_0.x() &&
           high_0.y() >= low_1.y() && high_1.y() >= low_0.y() &&
           high_0.z() >= low_1.z() && high_1.z() >= low_0.z();
}
