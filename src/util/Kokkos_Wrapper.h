//
// Created by alex on 8/26/24.
//

#ifndef KOMD_KOKKOS_WRAPPER_H
#define KOMD_KOKKOS_WRAPPER_H

#include "math/Array.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"
#include "Kokkos_ReductionIdentity.hpp"

namespace Kokkos {
    template<>
    struct reduction_identity<math::d3> {
        KOKKOS_FORCEINLINE_FUNCTION constexpr static math::d3 sum() { return {0, 0, 0}; }
        KOKKOS_FORCEINLINE_FUNCTION constexpr static math::d3 prod() { return {1, 1, 1}; }
        KOKKOS_FORCEINLINE_FUNCTION constexpr static math::d3 max() { return {-DBL_MAX, -DBL_MAX, -DBL_MAX}; }
        KOKKOS_FORCEINLINE_FUNCTION constexpr static math::d3 min() { return {DBL_MAX, DBL_MAX, DBL_MAX}; }
    };

    template<>
    KOKKOS_INLINE_FUNCTION void atomic_add<math::d3>(math::d3* const dest, desul::Impl::dont_deduce_this_parameter_t<const math::d3> val) {
        atomic_add(&dest->x(), val.x());
        atomic_add(&dest->y(), val.y());
        atomic_add(&dest->z(), val.z());
    }

    template<>
    KOKKOS_INLINE_FUNCTION void atomic_sub<math::d3>(math::d3* const dest, desul::Impl::dont_deduce_this_parameter_t<const math::d3> val) {
        atomic_sub(&dest->x(), val.x());
        atomic_sub(&dest->y(), val.y());
        atomic_sub(&dest->z(), val.z());
    }

    template<>
    KOKKOS_INLINE_FUNCTION void atomic_mul<math::d3>(math::d3* const dest, desul::Impl::dont_deduce_this_parameter_t<const math::d3> val) {
        atomic_mul(&dest->x(), val.x());
        atomic_mul(&dest->y(), val.y());
        atomic_mul(&dest->z(), val.z());
    }

    template<>
    KOKKOS_INLINE_FUNCTION void atomic_div<math::d3>(math::d3* const dest, desul::Impl::dont_deduce_this_parameter_t<const math::d3> val) {
        atomic_div(&dest->x(), val.x());
        atomic_div(&dest->y(), val.y());
        atomic_div(&dest->z(), val.z());
    }
};

/// Kokkos Wrapper
namespace KW {
#if defined(KOKKOS_ENABLE_CUDA)
    static constexpr bool enabled_cuda = true;
    using memory_space = Kokkos::SharedSpace;
    using contribution = Kokkos::Experimental::ScatterAtomic;
    using duplication = Kokkos::Experimental::ScatterNonDuplicated;
#else
    static constexpr bool enabled_cuda = false;
    using memory_space = Kokkos::HostSpace;
    using contribution = Kokkos::Experimental::ScatterNonAtomic;
    using duplication = Kokkos::Experimental::ScatterDuplicated;
#endif

    using exec_space = Kokkos::DefaultExecutionSpace;
    template<typename T>
    using vec_t = Kokkos::View<T*, memory_space>;

    template<typename T>
    using vec_scatter_t = Kokkos::Experimental::ScatterView<T*,
            exec_space::array_layout,
            exec_space,
            Kokkos::Experimental::ScatterSum,
            duplication,
            contribution>;
};

#endif //KOMD_KOKKOS_WRAPPER_H
