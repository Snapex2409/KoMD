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
m_use_soa(Registry::instance->configuration()->enableSOA),
m_enable_opt_approx(Registry::instance->configuration()->sensor_temp_opt),
m_temperature(0) { }

TemperatureSensor::TemperatureSensor(const math::d3 &low, const math::d3 &high) :
Sensor("Temperature"), p_low(low), p_high(high),
m_use_soa(Registry::instance->configuration()->enableSOA),
m_enable_opt_approx(Registry::instance->configuration()->sensor_temp_opt),
m_temperature(0) { }

void TemperatureSensor::measure() {
    auto config = Registry::instance->configuration();
    auto container = Registry::instance->moleculeContainer();
    auto& cells = container->getCells();
    const math::ul3 cell_dims = cells.dims();

    double mv2 = 0;
    double num_sites = 0;
    // loop over all non-halo cells
    for (uint64_t z = 1; z < cell_dims.z()-1; z++) {
        for (uint64_t y = 1; y < cell_dims.y()-1; y++) {
            for (uint64_t x = 1; x < cell_dims.x()-1; x++) {
                Cell& cell = cells[x, y, z];
                if (!math::boxesIntersect(cell.low(), cell.high(), p_low, p_high)) continue;

                if (!m_use_soa) {
                    for (Molecule& molecule : cell.molecules()) {
                        for (Site& site : molecule.getSites()) {
                            if (!math::pointInBox(site.r_arr(), p_low, p_high)) continue;
                            mv2 += site.getMass() * site.v_arr().dot(site.v_arr());
                            num_sites += 1;
                        }
                    }
                }
                else {
                    SOA& soa = cell.soa();
                    for (uint64_t idx = 0; idx < soa.size(); idx++) {
                        if (!math::pointInBox(soa.r()[idx], p_low, p_high)) continue;
                        mv2 += soa.mass()[idx] * soa.v()[idx].dot(soa.v()[idx]);
                        num_sites += 1;
                    }
                }
            }
        }
    }

    static const SciValue convDaInvkb = Constants::conv_Da_kg / Constants::kB;
    static const SciValue factor = Constants::conv_Aps_ms * Constants::conv_Aps_ms * convDaInvkb;
    m_temperature = factor * (mv2 / num_sites);
}

void TemperatureSensor::write(uint64_t simstep) { }

double TemperatureSensor::getTemperature() { return m_temperature; }
