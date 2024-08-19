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
        SOA& soa;
        const double stiffness_factor;
    };
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    double m_stiffness_factor;
};


#endif //KOMD_FENE_H
