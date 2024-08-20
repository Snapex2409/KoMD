//
// Created by alex on 8/8/24.
//

#ifndef KOMD_LJ12_6_SENSOR_H
#define KOMD_LJ12_6_SENSOR_H


#include "Potential_Sensor.h"
#include "container/SOA.h"
#include "math/Array.h"

/**
 * Measures potential of LJ126
 * */
class LJ12_6_Sensor : public Potential_Sensor {
public:
    LJ12_6_Sensor();
    virtual ~LJ12_6_Sensor() = default;
    void measure() override;

    /// @returns current total potential value
    double getCurrentPotential() { return m_total_pot[0]; }

    /// Kernel for device
    struct LJ12_6_Pot {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        /// shallow copy of SOA::id
        SOA::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r;
        /// shallow copy of SOA::sigma
        SOA::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        SOA::vec_t<double> eps;

        /// ScatterView for Potential_Sensor::p_data_u
        SOA::vec_scatter_t<double> data_u_scatter;
        /// ScatterView for Potential_Sensor::p_data_f
        SOA::vec_scatter_t<double> data_f_scatter;
        /// ScatterView for Potential_Sensor::p_count_u
        SOA::vec_scatter_t<uint64_t> count_u_scatter;
        /// ScatterView for Potential_Sensor::p_count_f
        SOA::vec_scatter_t<uint64_t> count_f_scatter;
        /// ScatterView for total potential buffer
        SOA::vec_scatter_t<double> total_pot_scatter;

        /// cutoff radius squared
        const double cutoff2;
        /// max sigma
        const double max_sigma;
        /// number of discretization bins
        const uint64_t bins;
    };

    /// Kernel for device
    struct LJ12_6_PotPair {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        /// shallow copy of SOA::id
        SOA::vec_t<uint64_t> id0;
        /// shallow copy of SOA::id
        SOA::vec_t<uint64_t> id1;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r0;
        /// shallow copy of SOA::r
        SOA::vec_t<math::d3> r1;
        /// shallow copy of SOA::sigma
        SOA::vec_t<double> sig0;
        /// shallow copy of SOA::sigma
        SOA::vec_t<double> sig1;
        /// shallow copy of SOA::epsilon
        SOA::vec_t<double> eps0;
        /// shallow copy of SOA::epsilon
        SOA::vec_t<double> eps1;

        /// ScatterView for Potential_Sensor::p_data_u
        SOA::vec_scatter_t<double> data_u_scatter;
        /// ScatterView for Potential_Sensor::p_data_f
        SOA::vec_scatter_t<double> data_f_scatter;
        /// ScatterView for Potential_Sensor::p_count_u
        SOA::vec_scatter_t<uint64_t> count_u_scatter;
        /// ScatterView for Potential_Sensor::p_count_f
        SOA::vec_scatter_t<uint64_t> count_f_scatter;
        /// ScatterView for total potential buffer
        SOA::vec_scatter_t<double> total_pot_scatter;

        /// cutoff radius squared
        const double cutoff2;
        /// max sigma
        const double max_sigma;
        /// number of discretization bins
        const uint64_t bins;
    };
protected:
    void handleCell(Cell &cell) override;
    void handleCellPair(Cell &cell0, Cell &cell1) override;
private:
    /// cutoff radius squared
    double m_cutoff2;
    /// total potential buffer
    SOA::vec_t<double> m_total_pot;
    /// scatter view for total potential buffer
    SOA::vec_scatter_t<double> m_total_pot_scatter;
};


#endif //KOMD_LJ12_6_SENSOR_H
