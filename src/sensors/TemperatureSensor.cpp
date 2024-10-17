//
// Created by alex on 8/8/24.
//

#include "TemperatureSensor.h"
#include "Registry.h"
#include "math/SciValue.h"
#include "math/Geometry.h"
#include "util/constants.h"

TemperatureSensor::TemperatureSensor() : Sensor("Temperature"),
p_low(Registry::instance->configuration()->domainLow),
p_high(Registry::instance->configuration()->domainHigh),
m_enable_opt_approx(Registry::instance->configuration()->sensor_temp_opt),
m_temperature(0) { }

TemperatureSensor::TemperatureSensor(const math::d3 &low, const math::d3 &high) :
Sensor("Temperature"), p_low(low), p_high(high),
m_enable_opt_approx(Registry::instance->configuration()->sensor_temp_opt),
m_temperature(0) { }

void TemperatureSensor::measure() {
    auto container = Registry::instance->moleculeContainer();
    SOA& soa = container->getSOA();

    auto& r = soa.r();
    auto& v = soa.v();
    auto& m = soa.mass();
    double num_sites = 0;
    double mv2 = 0;

    for (int idx = 0; idx < soa.size(); ++idx) {
        if (!math::pointInBox(r[idx], p_low, p_high)) continue;
        const auto vel = v[idx];
        mv2 += m[idx] * vel.dot(vel);
        num_sites += 1;
    }

    static const SciValue convDaInvkb = Constants::conv_Da_kg / Constants::kB;
    static const SciValue factor = Constants::conv_Aps_ms * Constants::conv_Aps_ms * convDaInvkb;
    num_sites *= 3; // for 3 dimensions
    m_temperature = factor * (mv2 / num_sites);
}

void TemperatureSensor::write(uint64_t simstep) { }

double TemperatureSensor::getTemperature() { return m_temperature; }

void TemperatureSensor::Temperature_Kernel::operator()(int idx) const {
    /*if (!math::pointInBox(r[idx], low, high)) return;
    const auto vel = v[idx];

    auto mv2_access = mv2.access();
    auto num_sites_access = num_sites.access();
    const auto e_kin_contrib = m[idx] * vel.dot(vel);
    if (e_kin_contrib != 0.0) mv2_access(0) += e_kin_contrib;
    num_sites_access(0) += 1;*/
}
