//
// Created by alex on 8/8/24.
//

#ifndef KOMD_POTENTIAL_SENSOR_H
#define KOMD_POTENTIAL_SENSOR_H


#include "Sensor.h"
#include "potentials/ForceFunctor.h"
#include "math/Array.h"

#include <vector>

class Potential_Sensor : public Sensor, protected ForceFunctor {
public:
    Potential_Sensor(std::string name, uint64_t bins);
    virtual ~Potential_Sensor() = default;

    void measure() override;

    void write(uint64_t simstep) override;

protected:
    uint64_t get_bin(double r);

    uint64_t p_bins;
    std::vector<double> p_data_u;
    std::vector<double> p_data_f;
    std::vector<uint64_t> p_count_u;
    std::vector<uint64_t> p_count_f;
private:
    const double m_sigma;
};


#endif //KOMD_POTENTIAL_SENSOR_H
