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
    struct Step0 {
        KOKKOS_FUNCTION void operator()(int idx) const;
        SOA::vec_t<math::d3> r;
        SOA::vec_t<math::d3> v;
        SOA::vec_t<math::d3> f;
        SOA::vec_t<double> m;
        const double dt_halve;
        const double dt;
    };

    struct Step1 {
        KOKKOS_FUNCTION void operator()(int idx) const;
        SOA::vec_t<math::d3> v;
        SOA::vec_t<math::d3> f;
        SOA::vec_t<double> m;
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
    double m_delta_t;
};


#endif //KOMD_INTEGRATOR_H
