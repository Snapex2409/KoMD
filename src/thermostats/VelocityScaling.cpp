//
// Created by alex on 8/10/24.
//

#include "VelocityScaling.h"
#include "Registry.h"
#include "sensors/TemperatureSensor.h"

VelocityScaling::VelocityScaling() : Thermostat(),
m_temp_target(Registry::instance->configuration()->temperature),
m_use_soa(Registry::instance->configuration()->enableSOA) { }

void VelocityScaling::apply() {
    auto temp_sensor = Registry::instance->temperature_sensor();
    const double temp_current = temp_sensor->getTemperature();
    const double beta = std::sqrt(m_temp_target / temp_current);

    auto config = Registry::instance->configuration();
    auto container = Registry::instance->moleculeContainer();
    auto& cells = container->getCells();
    const math::ul3 cell_dims = cells.dims();

    // loop over all non-halo cells
    for (uint64_t z = 1; z < cell_dims.z()-1; z++) {
        for (uint64_t y = 1; y < cell_dims.y()-1; y++) {
            for (uint64_t x = 1; x < cell_dims.x()-1; x++) {
                Cell& cell = cells[x, y, z];

                if (!m_use_soa) {
                    for (Molecule& molecule : cell.molecules()) {
                        for (Site& site : molecule.getSites()) {
                            site.v_arr() *= beta;
                        }
                    }
                }
                else {
                    SOA& soa = cell.soa();
                    for (uint64_t idx = 0; idx < soa.size(); idx++) {
                        soa.v()[idx] *= beta;
                    }
                }
            }
        }
    }
}

