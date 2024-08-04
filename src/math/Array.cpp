//
// Created by alex on 8/1/24.
//
#include "Array.h"

[[maybe_unused]] math::i3 math::ceil(const math::f3 &vec) {
    i3 result {0, 0, 0};
    result[0] = std::ceil(vec[0]);
    result[1] = std::ceil(vec[1]);
    result[2] = std::ceil(vec[2]);
    return result;
}

math::i3 math::floor(const math::f3 &vec) {
    i3 result {0, 0, 0};
    result[0] = std::floor(vec[0]);
    result[1] = std::floor(vec[1]);
    result[2] = std::floor(vec[2]);
    return result;
}

math::l3 math::ceil(const math::d3 &vec) {
    l3 result {0, 0, 0};
    result[0] = std::ceil(vec[0]);
    result[1] = std::ceil(vec[1]);
    result[2] = std::ceil(vec[2]);
    return result;
}

math::l3 math::floor(const math::d3 &vec) {
    l3 result {0, 0, 0};
    result[0] = std::floor(vec[0]);
    result[1] = std::floor(vec[1]);
    result[2] = std::floor(vec[2]);
    return result;
}

math::ui3 math::uceil(const math::f3 &vec) {
    ui3 result {0, 0, 0};
    result[0] = std::ceil(vec[0]);
    result[1] = std::ceil(vec[1]);
    result[2] = std::ceil(vec[2]);
    return result;
}

math::ui3 math::ufloor(const math::f3 &vec) {
    ui3 result {0, 0, 0};
    result[0] = std::floor(vec[0]);
    result[1] = std::floor(vec[1]);
    result[2] = std::floor(vec[2]);
    return result;
}

math::ul3 math::uceil(const math::d3 &vec) {
    ul3 result {0, 0, 0};
    result[0] = std::ceil(vec[0]);
    result[1] = std::ceil(vec[1]);
    result[2] = std::ceil(vec[2]);
    return result;
}

math::ul3 math::ufloor(const math::d3 &vec) {
    ul3 result {0, 0, 0};
    result[0] = std::floor(vec[0]);
    result[1] = std::floor(vec[1]);
    result[2] = std::floor(vec[2]);
    return result;
}