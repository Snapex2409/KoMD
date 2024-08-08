//
// Created by alex on 8/8/24.
//

#ifndef KOMD_SENSOR_H
#define KOMD_SENSOR_H

#include <cstdint>
#include <string>
#include <utility>

class Sensor {
public:
    explicit Sensor(std::string  name) : m_name(std::move(name)) {}
    virtual ~Sensor() = default;

    virtual void measure() = 0;
    virtual void write(uint64_t simstep) = 0;
    std::string name() { return m_name; }
private:
    std::string m_name;
};

#endif //KOMD_SENSOR_H
