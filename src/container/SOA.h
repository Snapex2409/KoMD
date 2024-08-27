//
// Created by alex on 7/31/24.
//

#ifndef KOMD_SOA_H
#define KOMD_SOA_H

#include <vector>
#include "math/Array.h"
#include "util/Kokkos_Wrapper.h"

class SOA {
public:
    explicit SOA();
    void resize(uint64_t newSize);
    void createBuffers(uint64_t size);

    KOKKOS_INLINE_FUNCTION KW::vec_t<double>& epsilon() { return m_epsilon; };
    KOKKOS_INLINE_FUNCTION KW::vec_t<double>& sigma() { return m_sigma; };
    KOKKOS_INLINE_FUNCTION KW::vec_t<double>& mass() { return m_mass; };
    KOKKOS_INLINE_FUNCTION KW::vec_t<uint64_t>& id() { return m_id; };
    KOKKOS_INLINE_FUNCTION KW::vec_t<uint64_t>& idx() { return m_idx; };
    KOKKOS_INLINE_FUNCTION KW::vec_t<math::d3>& r() { return m_r; };
    KOKKOS_INLINE_FUNCTION KW::vec_t<math::d3>& f() { return m_f; };
    KOKKOS_INLINE_FUNCTION KW::vec_t<math::d3>& v() { return m_v; };

    KOKKOS_INLINE_FUNCTION uint64_t size() { return m_current_size; }
private:
    KW::vec_t<double> m_epsilon;
    KW::vec_t<double> m_sigma;
    KW::vec_t<double> m_mass;
    KW::vec_t<uint64_t> m_id;
    KW::vec_t<uint64_t> m_idx;
    KW::vec_t<math::d3> m_r;
    KW::vec_t<math::d3> m_f;
    KW::vec_t<math::d3> m_v;

    uint64_t m_current_size;
};


#endif //KOMD_SOA_H
