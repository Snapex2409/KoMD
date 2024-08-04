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
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    inline void computeForce(Site& site0, Site& site1) const;
    inline void computeForceSOA(uint64_t idx_0, uint64_t idx_1,
                                SOA::vec_t<math::d3>& r0, SOA::vec_t<math::d3>& r1,
                                SOA::vec_t<math::d3>& f0, SOA::vec_t<math::d3>& f1,
                                SOA::vec_t<double>& sigmas0, SOA::vec_t<double>& sigmas1,
                                SOA::vec_t<double>& epsilons0, SOA::vec_t<double>& epsilons1) const;
    double m_cutoff2;
};


#endif //KOMD_LJ12_6_H
