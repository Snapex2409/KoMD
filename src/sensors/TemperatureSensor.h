//
// Created by alex on 8/8/24.
//

#ifndef KOMD_TEMPERATURE_SENSOR_H
#define KOMD_TEMPERATURE_SENSOR_H


#include "Sensor.h"
#include "math/Array.h"

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

protected:
    math::d3 p_low;
    math::d3 p_high;
private:
    bool m_enable_opt_approx;
    double m_temperature;
};

#endif //KOMD_TEMPERATURE_SENSOR_H
