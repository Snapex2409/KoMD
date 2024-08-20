//
// Created by alex on 8/4/24.
//

#ifndef KOMD_Limit_H
#define KOMD_Limit_H

#include "ForceFunctor.h"
#include "container/SOA.h"
#include <cstdint>

class Site;

class Limit : public ForceFunctor {
public:
    Limit();
    ~Limit() override = default;
    struct Limit_Force {
        KOKKOS_FUNCTION void operator()(int idx) const;
        SOA::vec_t<math::d3> f;
        SOA::vec_t<double> sig;
        SOA::vec_t<double> eps;
        const double limit_factor;
    };
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    double m_limit_factor;
};


#endif //KOMD_Limit_H
