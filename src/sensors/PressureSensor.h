//
// Created by alex on 8/8/24.
//

#ifndef KOMD_PRESSURE_H
#define KOMD_PRESSURE_H

#include "Sensor.h"
#include "util/Kokkos_Wrapper.h"
#include "potentials/ForceFunctor.h"
#include "potentials/ForceFunctor3B.h"

/**
* Pressure Measurement for single site particles
*/
class PressureSensor : public Sensor {
public:
    enum PressureDim {XX, XY, YY, XZ, YZ, ZZ, NUM_PRESSURES};
    PressureSensor();
    ~PressureSensor() override = default;

    void measure() override;

    void write(uint64_t simstep) override;

    /// momentum kernel
    struct MomentumKernel
    {
        KOKKOS_FUNCTION void operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy, double& acc_xz, double& acc_yz, double& acc_zz) const;
        /// particle velocities
        KW::vec_t<math::d3> v;
        /// particle masses
        KW::vec_t<double> m;
    };

    struct LJ12_6_Pressure {
        KOKKOS_FUNCTION void operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy, double& acc_xz, double& acc_yz, double& acc_zz) const;
        /// shallow copy of PairList::pairs
        KW::nvec_t<int, 2> pairs;
        /// shallow copy of PairList::pair_offsets
        KW::nvec_t<math::d3, 2> pair_offsets;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::id
        KW::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps;
        /// cutoff radius squared
        const double cutoff2;
    };

    struct ATM2B_Pressure {
        KOKKOS_FUNCTION void operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy, double& acc_xz, double& acc_yz, double& acc_zz) const;
        /// shallow copy of PairList::pairs
        KW::nvec_t<int, 2> pairs;
        /// shallow copy of PairList::pair_offsets
        KW::nvec_t<math::d3, 2> pair_offsets;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::id
        KW::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps;
        /// cutoff radius squared
        const double cutoff2;
        /// scaling factor
        const double potential_factor;
    };

    struct ATM_Pressure {
        KOKKOS_FUNCTION void operator()(int idx, double& acc_xx, double& acc_xy, double& acc_yy, double& acc_xz, double& acc_yz, double& acc_zz) const;
        /// shallow copy of TripleList::triplets
        KW::nvec_t<int, 3> triplets;
        /// shallow copy of TripleList::triple_offsets
        KW::nvec_t<math::d3, 3> triple_offsets;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// Energy parameter
        const double nu;
        /// cutoff radius squared
        const double cutoff2;
    };

    const KW::vec_t<double>& getPressureTensor() const { return m_pressures; }
private:
    void contributeMomentum();
    void contributePotentials();

    void runFF2(std::shared_ptr<ForceFunctor> ff);
    void runFF3(std::shared_ptr<ForceFunctor3B> ff);

    /// 3b nu param
    double m_nu;
    /// pressure matrix / stress tensor... idk what term is correct here
    KW::vec_t<double> m_pressures;
};


#endif //KOMD_PRESSURE_H
