//
// Created by alex on 7/31/24.
//

#include <iostream>
#include "ForceFunctor.h"
#include "Registry.h"

ForceFunctor::ForceFunctor() : p_run_cells(true), p_run_pairs(true), p_soa(Registry::instance->moleculeContainer()->getSOA()) {}

void ForceFunctor::operator()() {
    auto container = Registry::instance->moleculeContainer();

    // single cells
    if (p_run_cells) {
        // compute forces
        for (auto it = container->iteratorCell(); it->isValid(); ++(*it)) {
            Cell& cell = it->cell();
            handleCell(cell);
        }
        Kokkos::fence("Potential - Single Cells");
    }

    // cell pairs
    if (p_run_pairs) {
        // compute forces
        for (auto it = container->iteratorC08(); it->isValid(); ++(*it)) {
            if (it->colorSwitched()) Kokkos::fence("Potential - Cell Pairs Color");

            Cell& cell0 = it->cell0();
            Cell& cell1 = it->cell1();
            const math::d3 cell0_shift = it->getCell0Shift();
            const math::d3 cell1_shift = it->getCell1Shift();
            handleCellPair(cell0, cell1, cell0_shift, cell1_shift);
        }
    }
    Kokkos::fence("Potential - Cell Pairs End");
}
