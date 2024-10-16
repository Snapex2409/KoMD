//
// Created by alex on 10/15/24.
//

#ifndef ADR_FORCE_H
#define ADR_FORCE_H

#include "math/Array.h"
#include "potentials/ForceFunctor.h"

class ADR_Force : public ForceFunctor {
public:
    explicit ADR_Force(KW::vec_t<double>& weights);
    ~ADR_Force() override = default;

    /// Kernel for device
    struct Force_Kernel {
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
        /// shallow copy of SOA::mass
        KW::vec_t<double> mass;
        /// cutoff radius squared
        const double cutoff2;
        /// FENE stiffness
        const double stiffness_factor;
        /// ADR weights
        KW::vec_t<double> weights;

        KOKKOS_FUNCTION void handle_FENE(int idx, uint64_t s_idx_0, uint64_t s_idx_1) const;
        KOKKOS_FUNCTION void handle_LJ(int idx, uint64_t s_idx_0, uint64_t s_idx_1) const;
    };
protected:
    void handlePairList(PairList &pairList) override;

private:
    /// cutoff radius squared
    const double m_cutoff2;
    /// FENE stiffness factor
    const double m_stiffness_factor;
    /// ref to weights buffer
    KW::vec_t<double>& m_weights_ref;
};



#endif //ADR_FORCE_H
