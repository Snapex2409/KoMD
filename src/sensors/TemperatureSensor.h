//
// Created by alex on 8/8/24.
//

#ifndef KOMD_TEMPERATURE_SENSOR_H
#define KOMD_TEMPERATURE_SENSOR_H


#include "Sensor.h"
#include "util/Kokkos_Wrapper.h"

/**
 * Measures temperature using: sum mv2 = N*k_B*T \n
 * If opt approx is enabled, then temperature is measured following: https://pubs.acs.org/doi/10.1021/acs.jctc.8b00874
 * */
class TemperatureSensor : public Sensor {
public:
    /**
     * Measures in full domain
     * */
    TemperatureSensor();

    /**
     * Measures in selected region
     * */
    TemperatureSensor(const math::d3& low, const math::d3& high);

    ~TemperatureSensor() override = default;
    void measure() override;

    /// NOP
    void write(uint64_t simstep) override;

    /// @returns current measured temperature
    [[maybe_unused, nodiscard]] double getTemperature();

    /// Kernel for device
    struct Temperature_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        /// shallow copy of SOA::r
        KW::vec_t<math::d3> r;
        /// shallow copy of SOA::v
        KW::vec_t<math::d3> v;
        /// shallow copy of SOA::mass
        KW::vec_t<double> m;
        /// lower corner of measuring region
        const math::d3 low;
        /// upper corner of measuring region
        const math::d3 high;
        /// ScatterView to mv2 buffer
        KW::vec_scatter_t<double> mv2;
        /// ScatterView to num_sites buffer
        KW::vec_scatter_t<double> num_sites;
    };

protected:
    /// lower corner of measuring region
    math::d3 p_low;
    /// upper corner of measuring region
    math::d3 p_high;
private:
    bool m_enable_opt_approx;
    /// last measured temperature
    double m_temperature;
};

#endif //KOMD_TEMPERATURE_SENSOR_H
