//
// Created by alex on 8/4/24.
//

#ifndef KOMD_FENE_H
#define KOMD_FENE_H

#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"
#include <cstdint>

/**
 * Computes FENE between all sites of on molecule
 * */
class FENE : public ForceFunctor {
public:
    FENE();
    ~FENE() override = default;

    /// Kernel for device
    struct FENE_Force {
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
        /// shallow copy of Cell::indices
        KW::vec_t<uint64_t> indices;
        /// stiffness of FENE pot
        const double stiffness_factor;
    };
protected:
    void handleCell(Cell &cell) override;
    void handleCellPair(Cell &cell0, Cell &cell1, const math::d3& cell0_shift, const math::d3& cell1_shift) override;
private:
    /// stiffness of FENE pot
    double m_stiffness_factor;
};


#endif //KOMD_FENE_H
