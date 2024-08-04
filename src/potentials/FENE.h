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
    double m_stiffness_factor;
};


#endif //KOMD_FENE_H
