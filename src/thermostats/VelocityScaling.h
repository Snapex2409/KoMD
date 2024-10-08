//
// Created by alex on 8/10/24.
//

#ifndef VELOCITYSCALING_H
#define VELOCITYSCALING_H

#include "Thermostat.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Velocity scaling thermostat with no ramp
 * */
class VelocityScaling : public Thermostat {
public:
    VelocityScaling();
    ~VelocityScaling() override = default;

    /// scales the velocities of all non halo indices
    void apply() override;

    /// Kernel for device
    struct VS_Kernel {
        /**
         * Gets called from Kokkos
         * @param idx index of site inside current SOA
         * */
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of v vector
        KW::vec_t<math::d3> v;
        /// scaling factor
        const double beta;
    };
private:
    /// target temperature
    double m_temp_target;
};



#endif //VELOCITYSCALING_H
