//
// Created by alex on 8/4/24.
//

#ifndef KOMD_FENE_H
#define KOMD_FENE_H

#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Computes FENE between all sites of on molecule
 * */
class FENE : public ForceFunctor {
public:
    FENE();
    ~FENE() override = default;

    /// Kernel for device
    struct FENE_Force {
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
        /// stiffness of FENE pot
        const double stiffness_factor;
    };

protected:
    void handlePairList(PairList &pairList) override;

private:
    /// stiffness of FENE pot
    double m_stiffness_factor;
};


#endif //KOMD_FENE_H
