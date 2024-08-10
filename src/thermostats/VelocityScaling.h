//
// Created by alex on 8/10/24.
//

#ifndef VELOCITYSCALING_H
#define VELOCITYSCALING_H

#include "Thermostat.h"

class VelocityScaling : public Thermostat {
public:
    VelocityScaling();
    ~VelocityScaling() override = default;

    void apply() override;
private:
    double m_temp_target;
    bool m_use_soa;
};



#endif //VELOCITYSCALING_H
