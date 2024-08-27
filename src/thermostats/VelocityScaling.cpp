//
// Created by alex on 8/10/24.
//

#include "VelocityScaling.h"
#include "Registry.h"
#include "sensors/TemperatureSensor.h"

VelocityScaling::VelocityScaling() : Thermostat(),
m_temp_target(Registry::instance->configuration()->temperature) { }

void VelocityScaling::apply() {
    auto temp_sensor = Registry::instance->temperature_sensor();
    const double temp_current = temp_sensor->getTemperature();
    const double beta = std::sqrt(m_temp_target / temp_current);

    auto container = Registry::instance->moleculeContainer();
    SOA& soa = container->getSOA();
    Kokkos::parallel_for("Velocity Scaling", soa.size(), VS_Kernel(soa.v(), beta));
    Kokkos::fence("VS fence");
}

void VelocityScaling::VS_Kernel::operator()(int idx) const {
    v[idx] *= beta;
}
