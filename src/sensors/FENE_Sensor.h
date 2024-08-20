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

    struct FENE_Pot {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA::vec_t<uint64_t> id;
        SOA::vec_t<math::d3> r;
        SOA::vec_t<double> sig;
        SOA::vec_t<double> eps;
        SOA::vec_scatter_t<double> data_u_scatter;
        SOA::vec_scatter_t<double> data_f_scatter;
        SOA::vec_scatter_t<uint64_t> count_u_scatter;
        SOA::vec_scatter_t<uint64_t> count_f_scatter;
        const double stiffness_factor;
        const double max_sigma;
        const uint64_t bins;
    };
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    double m_stiffness_factor;
};


#endif //KOMD_FENE_SENSOR_H
