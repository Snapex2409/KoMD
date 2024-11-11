//
// Created by alex on 11/11/24.
//

#ifndef EXCLUISIONLJ_H
#define EXCLUISIONLJ_H

#include "potentials/ForceFunctor.h"

class ExcluisionLJ : public ForceFunctor {
public:
    ExcluisionLJ();
    ~ExcluisionLJ() override = default;

    struct ExLJ_Force {
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
        /// exclusion zone for reload: low
        math::d3 exclusion_low;
        /// exclusion zone for reload: high
        math::d3 exclusion_high;
    };

protected:
    void handlePairList(PairList &pairList) override;

private:
    /// cutoff radius squared
    double m_cutoff2;
    /// exclusion zone for reload: low
    math::d3 m_exclusion_low;
    /// exclusion zone for reload: high
    math::d3 m_exclusion_high;
};



#endif //EXCLUISIONLJ_H
