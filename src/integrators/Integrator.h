//
// Created by alex on 7/31/24.
//

#ifndef KOMD_INTEGRATOR_H
#define KOMD_INTEGRATOR_H

#include "container/SOA.h"
#include "util/Kokkos_Wrapper.h"

/**
 * 2 Step Leapfrog integration
 * */
class Integrator {
public:
    /// Kernel for device
    struct Step0 {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of SOA::v
        KW::vec_t<math::d3> v;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::mass
        KW::vec_t<double> m;
        /// delta_t / 2
        const double dt_halve;
        /// delta_t
        const double dt;
    };

    /// Kernel for device
    struct Step1 {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of SOA::v
        KW::vec_t<math::d3> v;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::mass
        KW::vec_t<double> m;
        /// delta_t / 2
        const double dt_halve;
    };

    /**
     * Creates new Integrator
     * */
    explicit Integrator();

    virtual ~Integrator() = default;

    /**
     * First half step of integration
     * */
    virtual void integrate0();

    /**
     * Second half step of integration
     * */
    virtual void integrate1();
protected:
    /// delta_t
    double p_delta_t;
};


#endif //KOMD_INTEGRATOR_H
