//
// Created by alex on 8/21/24.
//

#ifndef KOMD_OCLJ_H
#define KOMD_OCLJ_H

#include <cstdint>
#include "ForceFunctor.h"
#include "container/SOA.h"
#include "container/OneCell.h"

class OCLJ : public ForceFunctor {
public:
    OCLJ();
    ~OCLJ() override = default;

    /// Kernel for device
    struct OCLJ_Force {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1, int s_idx) const;
        /// ScatterView to SOA::f
        SOA::vec_scatter_t<math::d3> scatter;
        /// shallow copy of SOA::id
        SOA::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r;
        /// shallow copy of CoM
        SOA::vec_t<math::d3> com;
        /// shallow copy of SOA::sigma
        SOA::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        SOA::vec_t<double> eps;
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
