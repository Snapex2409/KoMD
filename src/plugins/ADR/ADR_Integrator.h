//
// Created by alex on 10/15/24.
//

#ifndef ADR_INTEGRATOR_H
#define ADR_INTEGRATOR_H

#include "integrators/Integrator.h"

class ADR_Integrator final : public Integrator{
public:
    explicit ADR_Integrator(double total_mass);

    ~ADR_Integrator() override = default;

    /// Kernel for device
    struct ADR_Step0 {
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
        /// CG replacement mass
        const double total_mass;
    };

    /// Kernel for device
    struct ADR_Step1 {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of SOA::v
        KW::vec_t<math::d3> v;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::mass
        KW::vec_t<double> m;
        /// delta_t / 2
        const double dt_halve;
        /// CG replacement mass
        const double total_mass;
    };

    void integrate0() override;

    void integrate1() override;
private:
    const double m_total_mass;
};



#endif //ADR_INTEGRATOR_H
