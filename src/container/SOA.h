//
// Created by alex on 7/31/24.
//

#ifndef KOMD_SOA_H
#define KOMD_SOA_H

#include <vector>
#include "math/Array.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

class SOA {
public:
    template<typename T>
    using vec_t = Kokkos::View<T*, Kokkos::SharedSpace>;
    template<typename T>
    using vec_scatter_t = Kokkos::Experimental::ScatterView<T*>;

    explicit SOA();
    void resize(uint64_t newSize);
    KOKKOS_INLINE_FUNCTION vec_scatter_t<double>& epsilonScatter() { return m_epsilon_scatter; };
    KOKKOS_INLINE_FUNCTION vec_scatter_t<double>& sigmaScatter() { return m_sigma_scatter; };
    KOKKOS_INLINE_FUNCTION vec_scatter_t<double>& massScatter() { return m_mass_scatter; };
    KOKKOS_INLINE_FUNCTION vec_scatter_t<uint64_t>& idScatter() { return m_id_scatter; };
    KOKKOS_INLINE_FUNCTION vec_scatter_t<math::d3>& rScatter() { return m_r_scatter; };
    KOKKOS_INLINE_FUNCTION vec_scatter_t<math::d3>& fScatter() { return m_f_scatter; };
    KOKKOS_INLINE_FUNCTION vec_scatter_t<math::d3>& vScatter() { return m_v_scatter; };

    KOKKOS_INLINE_FUNCTION vec_t<double>& epsilon() { return m_epsilon; };
    KOKKOS_INLINE_FUNCTION vec_t<double>& sigma() { return m_sigma; };
    KOKKOS_INLINE_FUNCTION vec_t<double>& mass() { return m_mass; };
    KOKKOS_INLINE_FUNCTION vec_t<uint64_t>& id() { return m_id; };
    KOKKOS_INLINE_FUNCTION vec_t<math::d3>& r() { return m_r; };
    KOKKOS_INLINE_FUNCTION vec_t<math::d3>& f() { return m_f; };
    KOKKOS_INLINE_FUNCTION vec_t<math::d3>& v() { return m_v; };

    KOKKOS_INLINE_FUNCTION uint64_t size() { return m_current_size; }
private:
    vec_scatter_t<double> m_epsilon_scatter;
    vec_scatter_t<double> m_sigma_scatter;
    vec_scatter_t<double> m_mass_scatter;
    vec_scatter_t<uint64_t> m_id_scatter;
    vec_scatter_t<math::d3> m_r_scatter;
    vec_scatter_t<math::d3> m_f_scatter;
    vec_scatter_t<math::d3> m_v_scatter;

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
