//
// Created by alex on 8/11/24.
//

#ifndef RDFSENSOR_H
#define RDFSENSOR_H

#include "Sensor.h"
#include <vector>

#include "Kokkos_Core.hpp"
#include "Kokkos_ScatterView.hpp"
#include "math/Array.h"

class Molecule;

class RDFSensor : public Sensor {
public:
    RDFSensor();
    virtual ~RDFSensor() = default;

    void measure() override;

    void write(uint64_t simstep) override;

    struct RDF_Kernel {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        Kokkos::View<math::d3*, Kokkos::SharedSpace>& positions;
        RDFSensor& sensor;
        const double max_r;
        const double delta_r;
        const uint64_t bins;
    };
private:
    static KOKKOS_INLINE_FUNCTION uint64_t getBin(double r, double delta_r, uint64_t bins) {
        if (r < 0) return 0;

        auto bin = static_cast<uint64_t>(r / delta_r);
        return Kokkos::clamp(bin, 0UL, bins - 1);
    }

    double m_max_r;
    double m_delta_r;
    double m_rho_0;
    Kokkos::View<double*, Kokkos::SharedSpace> m_bins;
    Kokkos::Experimental::ScatterView<double*> m_bins_scatter;
    double m_samples;
};



#endif //RDFSENSOR_H
