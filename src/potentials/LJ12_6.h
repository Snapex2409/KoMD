//
// Created by alex on 8/4/24.
//

#ifndef KOMD_LJ12_6_H
#define KOMD_LJ12_6_H

#include <cstdint>
#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Computes LJ126 between all sites of !different! molecules.
 * */
class LJ12_6 : public ForceFunctor {
public:
    LJ12_6();
    ~LJ12_6() override = default;

    /// Kernel for device
    struct LJ12_6_Force {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
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

    /// Kernel for device
    struct LJ12_6_ForcePair {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f0;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f1;

        /// shallow copy of SOA::id
        KW::vec_t<uint64_t> id0;
        /// shallow copy of SOA::id
        KW::vec_t<uint64_t> id1;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r0;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r1;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig0;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig1;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps0;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps1;
        /// cutoff radius squared
        const double cutoff2;
    };
protected:
    void handleCell(Cell &cell) override;
    void handleCellPair(Cell &cell0, Cell &cell1) override;

private:
    /// cutoff radius squared
    double m_cutoff2;
};


#endif //KOMD_LJ12_6_H
