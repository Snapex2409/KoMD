//
// Created by alex on 7/31/24.
//

#include "Integrator.h"
#include "Registry.h"
#include "molecule/Molecule.h"

Integrator::Integrator(double delta_t) :
m_delta_t(delta_t),
m_use_soa(Registry::instance->configuration()->enableSOA) { }

void Integrator::integrate0() {
    auto container = Registry::instance->moleculeContainer();
    auto& cells = container->getCells();
    const math::ul3 cell_dims = cells.dims();
    // loop over all non-halo cells
    const double dt_halve = m_delta_t * 0.5;
    for (uint64_t z = 1; z < cell_dims.z()-1; z++) {
        for (uint64_t y = 1; y < cell_dims.y()-1; y++) {
            for (uint64_t x = 1; x < cell_dims.x()-1; x++) {
                Cell& cell = cells[x, y, z];

                if (!m_use_soa) {
                    for (Molecule& molecule : cell.molecules()) {
                        for (Site& site : molecule.getSites()) {
                            const double dtInv2m = dt_halve / site.getMass();
                            site.v_arr() += site.f_arr() * dtInv2m;
                            site.r_arr() += site.v_arr() * m_delta_t;
                        }
                    }
                }
                else {
                    SOA& soa = cell.soa();
                    for (uint64_t idx = 0; idx < soa.size(); idx++) {
                        const double dtInv2m = dt_halve / soa.mass()[idx];
                        soa.v()[idx] += soa.f()[idx] * dtInv2m;
                        soa.r()[idx] += soa.v()[idx] * m_delta_t;
                    }
                }
            }
        }
    }
}

void Integrator::integrate1() {
    auto container = Registry::instance->moleculeContainer();
    auto& cells = container->getCells();
    const math::ul3 cell_dims = cells.dims();
    // loop over all non-halo cells
    const double dt_halve = m_delta_t * 0.5;
    for (uint64_t z = 1; z < cell_dims.z()-1; z++) {
        for (uint64_t y = 1; y < cell_dims.y()-1; y++) {
            for (uint64_t x = 1; x < cell_dims.x()-1; x++) {
                Cell& cell = cells[x, y, z];

                if (!m_use_soa) {
                    for (Molecule& molecule : cell.molecules()) {
                        for (Site& site : molecule.getSites()) {
                            const double dtInv2m = dt_halve / site.getMass();
                            site.v_arr() += site.f_arr() * dtInv2m;
                        }
                    }
                }
                else {
                    SOA& soa = cell.soa();
                    for (uint64_t idx = 0; idx < soa.size(); idx++) {
                        const double dtInv2m = dt_halve / soa.mass()[idx];
                        soa.v()[idx] += soa.f()[idx] * dtInv2m;
                    }
                }
            }
        }
    }
}
