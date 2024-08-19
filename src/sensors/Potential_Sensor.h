//
// Created by alex on 8/8/24.
//

#ifndef KOMD_POTENTIAL_SENSOR_H
#define KOMD_POTENTIAL_SENSOR_H


#include "Sensor.h"
#include "potentials/ForceFunctor.h"
#include "math/Array.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"

#include <vector>

class Potential_Sensor : public Sensor, protected ForceFunctor {
public:
    Potential_Sensor(std::string name, uint64_t bins);
    virtual ~Potential_Sensor() = default;

    void measure() override;

    void write(uint64_t simstep) override;

protected:
    static KOKKOS_INLINE_FUNCTION uint64_t get_bin(double r, double sigma, uint64_t bins) {
        if (r < 0) return 0;

        const double max_size = 3.0 * sigma;
        const double bin_width = max_size / static_cast<double>(bins);
        auto bin = static_cast<uint64_t>(r / bin_width);
        return Kokkos::clamp(bin, 0UL, bins-1);
    }

    uint64_t p_bins;
    Kokkos::View<double*, Kokkos::SharedSpace> p_data_u;
    Kokkos::View<double*, Kokkos::SharedSpace> p_data_f;
    Kokkos::View<uint64_t*, Kokkos::SharedSpace> p_count_u;
    Kokkos::View<uint64_t*, Kokkos::SharedSpace> p_count_f;

    Kokkos::Experimental::ScatterView<double*> p_data_u_scatter;
    Kokkos::Experimental::ScatterView<double*> p_data_f_scatter;
    Kokkos::Experimental::ScatterView<uint64_t*> p_count_u_scatter;
    Kokkos::Experimental::ScatterView<uint64_t*> p_count_f_scatter;
    const double p_sigma;
};


#endif //KOMD_POTENTIAL_SENSOR_H
