//
// Created by alex on 8/8/24.
//

#ifndef KOMD_SENSOR_H
#define KOMD_SENSOR_H

#include <cstdint>
#include <string>
#include <utility>

/**
 * Interface to measure any simulation property
 * */
class Sensor {
public:
    /**
     * Constructs a new Sensor
     * @param name of this sensor
     * */
    explicit Sensor(std::string  name) : m_name(std::move(name)) {}
    virtual ~Sensor() = default;

    /// measure the property
    virtual void measure() = 0;
    /// write data to file
    virtual void write(uint64_t simstep) = 0;
    /// @returns name of sensor
    std::string name() { return m_name; }
private:
    /// sensor name
    std::string m_name;
};

#endif //KOMD_SENSOR_H
