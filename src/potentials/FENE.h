//
// Created by alex on 8/4/24.
//

#ifndef KOMD_FENE_H
#define KOMD_FENE_H

#include "ForceFunctor.h"
#include "container/SOA.h"
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
        /// stiffness of FENE pot
        const double stiffness_factor;
    };
protected:
    void handleCell(Cell &cell) override;
    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    /// stiffness of FENE pot
    double m_stiffness_factor;
};


#endif //KOMD_FENE_H
