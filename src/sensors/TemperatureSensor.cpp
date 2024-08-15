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

    double mv2 = 0;
    double num_sites = 0;
    for (auto it = container->iterator(MoleculeContainer::SITE, MoleculeContainer::DOMAIN); it.isValid(); ++it) {
        if (!math::pointInBox(it.r(), p_low, p_high)) continue;
        mv2 += it.mass() * it.v().dot(it.v());
        num_sites += 1;
    }

    static const SciValue convDaInvkb = Constants::conv_Da_kg / Constants::kB;
    static const SciValue factor = Constants::conv_Aps_ms * Constants::conv_Aps_ms * convDaInvkb;
    num_sites *= 3; // for 3 dimensions
    m_temperature = factor * (mv2 / num_sites);
}

void TemperatureSensor::write(uint64_t simstep) { }

double TemperatureSensor::getTemperature() { return m_temperature; }
