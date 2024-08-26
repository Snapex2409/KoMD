//
// Created by alex on 8/21/24.
//

#ifndef KOMD_OCLJ_H
#define KOMD_OCLJ_H

#include <cstdint>
#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"
#include "container/OneCell.h"

class OCLJ : public ForceFunctor {
public:
    OCLJ();
    ~OCLJ() override = default;

    /// Kernel for device
    struct OCLJ_Force {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1, int s_idx) const;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::id
        KW::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of CoM
        KW::vec_t<math::d3> com;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps;
        /// cutoff radius squared
        const double cutoff2;
        /// domain size
        const math::d3 domain_size;
        /// lower threshold for shift
        const math::d3 low_bound;
    };
protected:
    void handleCell(Cell &cell) override;
    void handleCellPair(Cell&, Cell&) override { };

private:
    /// cutoff radius squared
    double m_cutoff2;
    /// domain size
    math::d3 m_domain_size;
    /// OC container
    std::shared_ptr<OneCell> m_container;
    /// lower threshold for shift
    math::d3 m_low_bound;
};


#endif //KOMD_OCLJ_H
