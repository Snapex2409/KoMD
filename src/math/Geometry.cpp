//
// Created by alex on 8/2/24.
//
#include "Geometry.h"

bool math::pointInBox(const math::d3 &point, const math::d3 &low, const math::d3 &high) {
    return point.x() >= low.x() && point.y() >= low.y() && point.z() >= low.z() &&
           point.x() <= high.x() && point.y() <= high.y() && point.z() <= high.z();
}
