//
// Created by alex on 8/22/24.
//

#ifndef KOMD_DISPLACEMENTSENSOR_H
#define KOMD_DISPLACEMENTSENSOR_H

#include "Sensor.h"
#include "container/SOA.h"

class DisplacementSensor : public Sensor {
public:
    DisplacementSensor();

    ~DisplacementSensor() override = default;

    void measure() override;

    void write(uint64_t simstep) override;

    double getDisplacement() { return m_lambda; }

    /// Kernel for device
    struct Displacement_Kernel {
        KOKKOS_FUNCTION void operator()(int idx, double& lx, double& ly, double& lz) const;
        /// center of mass positions
        SOA::vec_t<math::d3> com;
        /// size of a single lattice cell
        math::d3 lattice_unit_size;
    };
private:
    /// displacement factor
    double m_lambda;
    /// size of a single lattice cell
    math::d3 m_lattice_unit_size;
    /// com buffer
    SOA::vec_t<math::d3> m_com;
};


#endif //KOMD_DISPLACEMENTSENSOR_H
