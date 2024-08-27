//
// Created by alex on 8/4/24.
//

#ifndef KOMD_Limit_H
#define KOMD_Limit_H

#include "ForceFunctor.h"
#include "util/Kokkos_Wrapper.h"
#include <cstdint>

/**
 * Limits the values in force buffers to a specified absolute threshold. Is applied per site.
 * */
class Limit : public ForceFunctor {
public:
    Limit();
    ~Limit() override = default;

    /// Kernel for device
    struct Limit_Force {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of SOA::f
        KW::vec_t<math::d3> f;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps;
        /// shallow copy of Cell::indices
        KW::vec_t<uint64_t> indices;
        /// threshold
        const double limit_factor;
    };
protected:
    void handleCell(Cell &cell) override;
    void handleCellPair(Cell &cell0, Cell &cell1, const math::d3& cell1_shift) override;

private:
    /// threshold
    double m_limit_factor;
};


#endif //KOMD_Limit_H
