//
// Created by alex on 7/31/24.
//

#ifndef KOMD_INTEGRATOR_H
#define KOMD_INTEGRATOR_H

#include "container/SOA.h"
#include "Kokkos_Core.hpp"

/**
 * 2 Step Leapfrog integration
 * */
class Integrator {
public:
    /// Kernel for device
    struct Step0 {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r;
        /// shallow copy of SOA::v
        SOA::vec_t<math::d3> v;
        /// shallow copy of SOA::f
        SOA::vec_t<math::d3> f;
        /// shallow copy of SOA::mass
        SOA::vec_t<double> m;
        /// delta_t / 2
        const double dt_halve;
        /// delta_t
        const double dt;
    };

    /// Kernel for device
    struct Step1 {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of SOA::v
        SOA::vec_t<math::d3> v;
        /// shallow copy of SOA::f
        SOA::vec_t<math::d3> f;
        /// shallow copy of SOA::mass
        SOA::vec_t<double> m;
        /// delta_t / 2
        const double dt_halve;
    };

    /**
     * Creates new Integrator
     * */
    explicit Integrator();

    /**
     * First half step of integration
     * */
    void integrate0();

    /**
     * Second half step of integration
     * */
    void integrate1();
private:
    /// delta_t
    double m_delta_t;
};


#endif //KOMD_INTEGRATOR_H
