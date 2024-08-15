//
// Created by alex on 7/31/24.
//

#ifndef KOMD_SOA_H
#define KOMD_SOA_H

#include <vector>
#include "math/Array.h"

#include "Kokkos_Core.hpp"

class SOA {
public:
    template<typename T>
    using vec_t = Kokkos::View<T*, Kokkos::SharedSpace>;

    explicit SOA(uint64_t size = 0);
    void resize(uint64_t newSize);
    vec_t<double>& epsilon() { return m_epsilon; };
    vec_t<double>& sigma() { return m_sigma; };
    vec_t<double>& mass() { return m_mass; };
    vec_t<uint64_t>& id() { return m_id; };
    vec_t<math::d3>& r() { return m_r; };
    vec_t<math::d3>& f() { return m_f; };
    vec_t<math::d3>& v() { return m_v; };
    uint64_t size() { return m_current_size; }
private:
    vec_t<double> m_epsilon;
    vec_t<double> m_sigma;
    vec_t<double> m_mass;
    vec_t<uint64_t> m_id;
    vec_t<math::d3> m_r;
    vec_t<math::d3> m_f;
    vec_t<math::d3> m_v;

    uint64_t m_current_size;
};


#endif //KOMD_SOA_H
