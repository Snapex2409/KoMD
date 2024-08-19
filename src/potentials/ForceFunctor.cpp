//
// Created by alex on 7/31/24.
//

#include "ForceFunctor.h"
#include "Registry.h"

#include "Kokkos_Core.hpp"

ForceFunctor::ForceFunctor() : p_run_cells(true), p_run_pairs(true), p_run_contribution(true) {}

void ForceFunctor::operator()() {
    auto container = Registry::instance->moleculeContainer();

    // single cells
    if (p_run_cells) {
        for (auto it = container->iteratorCell(MoleculeContainer::DOMAIN); it.isValid(); ++it) {
            Cell& cell = it.cell();
            handleCell(cell);
            // was scatter view used?
            if (p_run_contribution) Kokkos::Experimental::contribute(cell.soa().f(), cell.soa().fScatter());
        }
        Kokkos::fence("Potential - Single Cells");
    }

    // cell pairs
    if (p_run_pairs) {
        for (auto it = container->iteratorC08(); it.isValid(); ++it) {
            if (it.colorSwitched()) Kokkos::fence("Potential - Cell Pairs Color");

            Cell& cell0 = it.cell0();
            Cell& cell1 = it.cell1();
            handleCellPair(cell0, cell1);
        }

        // was scatter view used?
        if (p_run_contribution) {
            Kokkos::fence("Potential - Cell Pairs Computation");
            for (auto it = container->iteratorCell(MoleculeContainer::DOMAIN); it.isValid(); ++it) {
                Cell& cell = it.cell();
                Kokkos::Experimental::contribute(cell.soa().f(), cell.soa().fScatter());
            }
        }
    }
    Kokkos::fence("Potential - Cell Pairs End");
}
