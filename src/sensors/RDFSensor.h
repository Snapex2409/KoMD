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
#include "container/SOA.h"

/**
 * Measures RDF
 * */
class RDFSensor : public Sensor {
public:
    RDFSensor();
    virtual ~RDFSensor() = default;
    void measure() override;
    void write(uint64_t simstep) override;

    /// Kernel for device
    struct RDF_Kernel {
        KOKKOS_FUNCTION void operator()(int idx_0, int idx_1) const;
        /// shallow copy of CoM positions
        SOA::vec_t<math::d3> positions;
        /// shallow copy of rdf bins
        SOA::vec_scatter_t<double> bins_scatter;
        /// maximum measure distance
        const double max_r;
        /// bin size
        const double delta_r;
        /// number of total bins
        const uint64_t bins;
    };
private:
    /**
     * @param r distance r in U(r)
     * @param delta_r size of one bin
     * @param bins number of bins
     * */
    static KOKKOS_INLINE_FUNCTION uint64_t getBin(double r, double delta_r, uint64_t bins) {
        if (r < 0) return bins-1;

        auto bin = static_cast<uint64_t>(r / delta_r);
        return Kokkos::clamp(bin, 0UL, bins - 1);
    }

    /// maximum measure distance
    double m_max_r;
    /// bin size
    double m_delta_r;
    /// bulk density
    double m_rho_0;
    /// rdf buffer
    SOA::vec_t<double> m_bins;
    /// ScatterView for rdf buffer
    SOA::vec_scatter_t<double> m_bins_scatter;
    /// number of samples
    double m_samples;
    /// com buffer
    SOA::vec_t<math::d3> m_com;
};



#endif //RDFSENSOR_H
