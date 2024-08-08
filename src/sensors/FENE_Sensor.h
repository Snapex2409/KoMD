//
// Created by alex on 8/8/24.
//

#ifndef KOMD_FENE_SENSOR_H
#define KOMD_FENE_SENSOR_H

#include "Potential_Sensor.h"
#include "container/SOA.h"
#include "math/Array.h"

class Site;

class FENE_Sensor : public Potential_Sensor {
public:
    FENE_Sensor();
    virtual ~FENE_Sensor() = default;

protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    inline void computeForce(Site& site0, Site& site1);
    inline void computeForceSOA(uint64_t idx_0, uint64_t idx_1,
                                SOA::vec_t<math::d3>& r0, SOA::vec_t<math::d3>& r1,
                                SOA::vec_t<math::d3>& f0, SOA::vec_t<math::d3>& f1,
                                SOA::vec_t<double>& sigmas0, SOA::vec_t<double>& sigmas1,
                                SOA::vec_t<double>& epsilons0, SOA::vec_t<double>& epsilons1);
    double m_stiffness_factor;
};


#endif //KOMD_FENE_SENSOR_H
