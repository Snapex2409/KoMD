//
// Created by alex on 8/8/24.
//

#ifndef KOMD_TEMPERATURE_SENSOR_H
#define KOMD_TEMPERATURE_SENSOR_H


#include "Sensor.h"
#include "math/Array.h"
#include "Kokkos_ScatterView.hpp"
#include "container/SOA.h"

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

    [[maybe_unused, nodiscard]] double getTemperature();

    struct Temperature_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        SOA::vec_t<math::d3> r;
        SOA::vec_t<math::d3> v;
        SOA::vec_t<double> m;
        const math::d3 low;
        const math::d3 high;
        SOA::vec_scatter_t<double> mv2;
        SOA::vec_scatter_t<double> num_sites;
    };
protected:
    math::d3 p_low;
    math::d3 p_high;
private:
    bool m_enable_opt_approx;
    double m_temperature;
};

#endif //KOMD_TEMPERATURE_SENSOR_H
