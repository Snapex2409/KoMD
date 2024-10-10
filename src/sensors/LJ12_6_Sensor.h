//
// Created by alex on 8/8/24.
//

#ifndef KOMD_LJ12_6_SENSOR_H
#define KOMD_LJ12_6_SENSOR_H


#include "Potential_Sensor.h"
#include "util/Kokkos_Wrapper.h"
#include <vector>

/**
 * Measures potential of LJ126
 * */
class LJ12_6_Sensor : public Potential_Sensor {
public:
    LJ12_6_Sensor();
    virtual ~LJ12_6_Sensor() = default;
    void measure() override;
    void write(uint64_t simstep) override;

    /// @returns current total potential value
    double getCurrentPotential() { return m_total_pot[0]; }

    /// Kernel for device
    struct LJ12_6_Pot {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of PairList::pairs
        KW::nvec_t<int, 2> pairs;
        /// shallow copy of PairList::pair_offsets
        KW::nvec_t<math::d3, 2> pair_offsets;
        /// shallow copy of SOA::id
        KW::vec_t<uint64_t> id;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of SOA::sigma
        KW::vec_t<double> sig;
        /// shallow copy of SOA::epsilon
        KW::vec_t<double> eps;

        /// ScatterView for Potential_Sensor::p_data_u
        KW::vec_scatter_t<double> data_u_scatter;
        /// ScatterView for Potential_Sensor::p_data_f
        KW::vec_scatter_t<double> data_f_scatter;
        /// ScatterView for Potential_Sensor::p_count_u
        KW::vec_scatter_t<uint64_t> count_u_scatter;
        /// ScatterView for Potential_Sensor::p_count_f
        KW::vec_scatter_t<uint64_t> count_f_scatter;
        /// ScatterView for total potential buffer
        KW::vec_scatter_t<double> total_pot_scatter;

        /// cutoff radius squared
        const double cutoff2;
        /// max sigma
        const double max_sigma;
        /// number of discretization bins
        const uint64_t bins;
    };

protected:
    void handlePairList(PairList &pairList) override;
private:
    /// cutoff radius squared
    double m_cutoff2;
    /// total potential buffer
    KW::vec_t<double> m_total_pot;
    /// scatter view for total potential buffer
    KW::vec_scatter_t<double> m_total_pot_scatter;
    /// history of global potential values
    std::vector<double> m_pot_hist;
};


#endif //KOMD_LJ12_6_SENSOR_H
