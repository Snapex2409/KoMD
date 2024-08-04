//
// Created by alex on 8/1/24.
//

#ifndef KOMD_VEC3D_H
#define KOMD_VEC3D_H


#include <vector>
#include "math/Array.h"

template<typename T>
class Vec3D {
public:
    /**
     * Constructs a 3D vector with no memory allocated
     * */
    Vec3D() = default;

    /**
     * Initializes to 3D vector with the provided dimensions.
     * */
    void init(const math::ul3& dim) {
        m_dim = dim;
        m_data.resize(dim.product());
    }

    /**
     * Constructs a 3D vector with the provided dimensions.
     * */
    explicit Vec3D(const math::ul3& dim) : m_dim(dim) {
        m_data.resize(dim.product());
    }

    /**
     * Access to data using coord tuple
     * */
    T& operator[](const math::ul3& coord) {
        return m_data[mapSpatial2Flat(coord)];
    }

    /**
     * Access to data using coord tuple
     * */
    T& operator[](uint64_t x, uint64_t y, uint64_t z) {
        return m_data[mapSpatial2Flat(x, y, z)];
    }

    /**
     * Access to data using coord tuple
     * */
    T operator[](const math::ul3& coord) const {
        return m_data[mapSpatial2Flat(coord)];
    }

    /**
     * Access to data using coord tuple
     * */
    T operator[](uint64_t x, uint64_t y, uint64_t z) const {
        return m_data[mapSpatial2Flat(x, y, z)];
    }

    /**
     * Access to raw pointer
     * */
    T* data() { return m_data.data(); }

    /**
     * Returns cref to dimensions
     * */
    [[nodiscard]] const math::ul3& dims() const { return m_dim; }
private:
    inline uint64_t mapSpatial2Flat(const math::ul3& coord) {
        return coord.x() + m_dim.x() * coord.y() + m_dim.x() * m_dim.y() * coord.z();
    }
    inline uint64_t mapSpatial2Flat(uint64_t x, uint64_t y, uint64_t z) {
        return x + m_dim.x() * y + m_dim.x() * m_dim.y() * z;
    }

    std::vector<T> m_data;
    math::ul3 m_dim;
};


#endif //KOMD_VEC3D_H
