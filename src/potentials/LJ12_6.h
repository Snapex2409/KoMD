//
// Created by alex on 8/4/24.
//

#ifndef KOMD_LJ12_6_H
#define KOMD_LJ12_6_H

#include <cstdint>
#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Computes LJ126 between all sites of !different! indices.
 * */
class LJ12_6 : public ForceFunctor {
public:
    LJ12_6();
    ~LJ12_6() override = default;

    /// Kernel for device
    struct LJ12_6_Force {
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
    };

protected:
    void handlePairList(PairList &pairList) override;

private:
    /// cutoff radius squared
    double m_cutoff2;
};


#endif //KOMD_LJ12_6_H
