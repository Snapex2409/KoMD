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

        /// stiffness of FENE pot
        const double stiffness_factor;
        /// max sigma
        const double max_sigma;
        /// number of discretization bins
        const uint64_t bins;
    };
protected:
    void handlePairList(PairList &pairList) override;
private:
    /// stiffness of FENE pot
    double m_stiffness_factor;
};


#endif //KOMD_FENE_SENSOR_H
