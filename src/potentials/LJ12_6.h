//
// Created by alex on 8/4/24.
//

#ifndef KOMD_LJ12_6_H
#define KOMD_LJ12_6_H

#include <cstdint>
#include "ForceFunctor.h"
#include "container/SOA.h"

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
        /// ScatterView to SOA::f
        SOA::vec_scatter_t<math::d3> scatter;
        /// shallow copy of SOA::id
        SOA::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r;
        /// shallow copy of SOA::sigma
        SOA::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        SOA::vec_t<double> eps;
        /// cutoff radius squared
        const double cutoff2;
    };

    /// Kernel for device
    struct LJ12_6_ForcePair {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        /// ScatterView to SOA::f
        SOA::vec_scatter_t<math::d3> scatter0;
        /// ScatterView to SOA::f
        SOA::vec_scatter_t<math::d3> scatter1;

        /// shallow copy of SOA::id
        SOA::vec_t<uint64_t> id0;
        /// shallow copy of SOA::id
        SOA::vec_t<uint64_t> id1;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r0;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r1;
        /// shallow copy of SOA::sigma
        SOA::vec_t<double> sig0;
        /// shallow copy of SOA::sigma
        SOA::vec_t<double> sig1;
        /// shallow copy of SOA::epsilon
        SOA::vec_t<double> eps0;
        /// shallow copy of SOA::epsilon
        SOA::vec_t<double> eps1;
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
