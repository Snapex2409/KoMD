//
// Created by alex on 8/10/24.
//

#ifndef VELOCITYSCALING_H
#define VELOCITYSCALING_H

#include "Thermostat.h"
#include "container/SOA.h"
#include "Kokkos_Core.hpp"


class VelocityScaling : public Thermostat {
public:
    VelocityScaling();
    ~VelocityScaling() override = default;

    void apply() override;

    struct VS_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        SOA::vec_t<math::d3> v;
        const double beta;
    };
private:
    double m_temp_target;
};



#endif //VELOCITYSCALING_H
