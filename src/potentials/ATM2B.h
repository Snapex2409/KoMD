//
// Created by alex on 8/4/24.
//

#ifndef KOMD_ATM2B_H
#define KOMD_ATM2B_H

#include <cstdint>
#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Computes LJ126 between all sites of !different! indices.
 * */
class ATM2B : public ForceFunctor {
public:
    ATM2B();
    ~ATM2B() override = default;

    /// Kernel for device
    struct ATM2B_Force {
        KOKKOS_FUNCTION void operator()(int idx) const;
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

    double getPotentialFactor() const { return m_potential_factor; }

protected:
    void handlePairList(PairList &pairList) override;

private:
    /// cutoff radius squared
    double m_cutoff2;
    /// scaling factor
    double m_potential_factor;
};


#endif //KOMD_ATM2B_H
