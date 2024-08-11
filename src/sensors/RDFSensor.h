//
// Created by alex on 8/11/24.
//

#ifndef RDFSENSOR_H
#define RDFSENSOR_H

#include "Sensor.h"
#include <vector>

class RDFSensor : public Sensor {
public:
    RDFSensor();
    virtual ~RDFSensor() = default;

    void measure() override;

    void write(uint64_t simstep) override;

private:
    uint64_t getBin(double r) const;
    double m_max_r;
    double m_delta_r;
    double m_rho_0;
    std::vector<double> m_bins;
    double m_samples;
};



#endif //RDFSENSOR_H
