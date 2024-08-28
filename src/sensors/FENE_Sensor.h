//
// Created by alex on 8/8/24.
//

#ifndef KOMD_FENE_SENSOR_H
#define KOMD_FENE_SENSOR_H

#include "Potential_Sensor.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Measures potential of FENE
 * */
class FENE_Sensor : public Potential_Sensor {
public:
    FENE_Sensor();
    virtual ~FENE_Sensor() = default;

    /// Kernel for device
    struct FENE_Pot {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        /// shallow copy of SOA::id
        KW::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps;
        /// shallow copy of Cell::indices
        KW::vec_t<uint64_t> indices;

        /// ScatterView for Potential_Sensor::p_data_u
        KW::vec_scatter_t<double> data_u_scatter;
        /// ScatterView for Potential_Sensor::p_data_f
        KW::vec_scatter_t<double> data_f_scatter;
        /// ScatterView for Potential_Sensor::p_count_u
        KW::vec_scatter_t<uint64_t> count_u_scatter;
        /// ScatterView for Potential_Sensor::p_count_f
        KW::vec_scatter_t<uint64_t> count_f_scatter;

        /// stiffness of FENE pot
        const double stiffness_factor;
        /// max sigma
        const double max_sigma;
        /// number of discretization bins
        const uint64_t bins;
    };
protected:
    void handleCell(Cell &cell) override;
    void handleCellPair(Cell &cell0, Cell &cell1, const math::d3& cell0_shift, const math::d3& cell1_shift) override;
private:
    /// stiffness of FENE pot
    double m_stiffness_factor;
};


#endif //KOMD_FENE_SENSOR_H
