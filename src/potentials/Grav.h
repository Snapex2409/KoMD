//
// Created by alex on 8/4/24.
//

#ifndef KOMD_GRAV_H
#define KOMD_GRAV_H

#include <cstdint>
#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Computes gravity between all sites of !different! indices.
 * */
class Grav : public ForceFunctor {
public:
    Grav();
    ~Grav() override = default;

    /// Kernel for device
    struct Grav_Force {
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
        /// shallow copy of SOA::mass
        KW::vec_t<double> m;
    };

protected:
    void handlePairList(PairList &pairList) override;
};


#endif //KOMD_GRAV_H
