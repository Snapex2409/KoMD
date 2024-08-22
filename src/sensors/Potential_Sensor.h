//
// Created by alex on 8/8/24.
//

#ifndef KOMD_POTENTIAL_SENSOR_H
#define KOMD_POTENTIAL_SENSOR_H


#include "Sensor.h"
#include "potentials/ForceFunctor.h"
#include "math/Array.h"
#include "container/SOA.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

#include <vector>

/**
 * Sensor to measure potentials
 * */
class Potential_Sensor : public Sensor, protected ForceFunctor {
public:
    /**
     * @param name of potential
     * @param bins number of bins for spatial discretization of U(r)
     * */
    Potential_Sensor(std::string name, uint64_t bins);
    virtual ~Potential_Sensor() = default;

    /// measure potential, force and write to buffer
    void measure() override;

    /// write to file
    void write(uint64_t simstep) override;

protected:
    /**
     * @param r distance r in U(r)
     * @param sigma maximum sigma value
     * @param bins number of bins for spatial discretization
     * */
    static KOKKOS_INLINE_FUNCTION uint64_t get_bin(double r, double sigma, uint64_t bins) {
        if (r < 0) return 0;

        const double max_size = 3.0 * sigma;
        const double bin_width = max_size / static_cast<double>(bins);
        auto bin = static_cast<uint64_t>(r / bin_width);
        return Kokkos::clamp(bin, 0UL, bins-1);
    }

    /// number of bins for spatial discretization of U(r)
    uint64_t p_bins;
    /// potential bins
    SOA::vec_t<double> p_data_u;
    /// force bins
    SOA::vec_t<double> p_data_f;
    /// potential number samples bins
    SOA::vec_t<uint64_t> p_count_u;
    /// force number samples bins
    SOA::vec_t<uint64_t> p_count_f;

    /// ScatterView for potential bins
    SOA::vec_scatter_t<double> p_data_u_scatter;
    /// ScatterView for force bins
    SOA::vec_scatter_t<double> p_data_f_scatter;
    /// ScatterView for potential number samples bins
    SOA::vec_scatter_t<uint64_t> p_count_u_scatter;
    /// ScatterView for force number samples bins
    SOA::vec_scatter_t<uint64_t> p_count_f_scatter;
    /// max sigma
    const double p_max_sigma;
};


#endif //KOMD_POTENTIAL_SENSOR_H
