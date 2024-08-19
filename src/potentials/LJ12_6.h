//
// Created by alex on 8/4/24.
//

#ifndef KOMD_LJ12_6_H
#define KOMD_LJ12_6_H

#include <cstdint>
#include "ForceFunctor.h"
#include "container/SOA.h"

class Site;

class LJ12_6 : public ForceFunctor {
public:
    LJ12_6();
    ~LJ12_6() override = default;

    struct LJ12_6_Force {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA& soa;
        const double cutoff2;
    };

    struct LJ12_6_ForcePair {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA& soa0;
        SOA& soa1;
        const double cutoff2;
    };
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    double m_cutoff2;
};


#endif //KOMD_LJ12_6_H
