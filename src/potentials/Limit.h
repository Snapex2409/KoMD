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
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    inline void computeForce(Site& site) const;
    inline void computeForceSOA(uint64_t idx, SOA::vec_t<math::d3>& r0, SOA::vec_t<math::d3>& f0, SOA::vec_t<double>& sigmas, SOA::vec_t<double>& epsilons) const;
    double m_limit_factor;
};


#endif //KOMD_Limit_H
