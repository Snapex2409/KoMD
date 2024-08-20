//
// Created by alex on 8/4/24.
//

#ifndef KOMD_FENE_H
#define KOMD_FENE_H

#include "ForceFunctor.h"
#include "container/SOA.h"
#include <cstdint>

class Site;

class FENE : public ForceFunctor {
public:
    FENE();
    ~FENE() override = default;
    struct FENE_Force {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA::vec_scatter_t<math::d3> scatter;
        SOA::vec_t<uint64_t> id;
        SOA::vec_t<math::d3> r;
        SOA::vec_t<double> sig;
        SOA::vec_t<double> eps;
        const double stiffness_factor;
    };
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    double m_stiffness_factor;
};


#endif //KOMD_FENE_H
