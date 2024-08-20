//
// Created by alex on 8/8/24.
//

#ifndef KOMD_LJ12_6_SENSOR_H
#define KOMD_LJ12_6_SENSOR_H


#include "Potential_Sensor.h"
#include "container/SOA.h"
#include "math/Array.h"

class Site;


class LJ12_6_Sensor : public Potential_Sensor {
public:
    LJ12_6_Sensor();
    virtual ~LJ12_6_Sensor() = default;
    void measure() override;
    double getCurrentPotential() { return m_total_pot[0]; }

    struct LJ12_6_Pot {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA::vec_t<uint64_t> id;
        SOA::vec_t<math::d3> r;
        SOA::vec_t<double> sig;
        SOA::vec_t<double> eps;

        SOA::vec_scatter_t<double> data_u_scatter;
        SOA::vec_scatter_t<double> data_f_scatter;
        SOA::vec_scatter_t<uint64_t> count_u_scatter;
        SOA::vec_scatter_t<uint64_t> count_f_scatter;
        SOA::vec_scatter_t<double> total_pot_scatter;
        const double cutoff2;
        const double max_sigma;
        const uint64_t bins;
    };

    struct LJ12_6_PotPair {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        SOA::vec_t<uint64_t> id0;
        SOA::vec_t<uint64_t> id1;
        SOA::vec_t<math::d3> r0;
        SOA::vec_t<math::d3> r1;
        SOA::vec_t<double> sig0;
        SOA::vec_t<double> sig1;
        SOA::vec_t<double> eps0;
        SOA::vec_t<double> eps1;

        SOA::vec_scatter_t<double> data_u_scatter;
        SOA::vec_scatter_t<double> data_f_scatter;
        SOA::vec_scatter_t<uint64_t> count_u_scatter;
        SOA::vec_scatter_t<uint64_t> count_f_scatter;
        SOA::vec_scatter_t<double> total_pot_scatter;
        const double cutoff2;
        const double max_sigma;
        const uint64_t bins;
    };
protected:
    void handleCell(Cell &cell) override;

    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    double m_cutoff2;
    Kokkos::View<double*, Kokkos::SharedSpace> m_total_pot;
    Kokkos::Experimental::ScatterView<double*> m_total_pot_scatter;
};


#endif //KOMD_LJ12_6_SENSOR_H
