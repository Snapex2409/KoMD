//
// Created by alex on 7/31/24.
//

#ifndef KOMD_FORCEFUNCTOR_H
#define KOMD_FORCEFUNCTOR_H

class Cell;
class SOA;
#include "math/Array.h"

/**
 * General force function interface
 * */
class ForceFunctor {
public:
    ForceFunctor();
    virtual ~ForceFunctor() = default;

    /**
     * This should compute some sort of forces.
     * */
    void operator()();

protected:
    /**
     * Handles a single cell
     * */
    virtual void handleCell(Cell& cell) = 0;

    /**
     * Handles a cell pair
     * */
    virtual void handleCellPair(Cell& cell0, Cell& cell1, const math::d3& cell0_shift, const math::d3& cell1_shift) = 0;

    /// Flag if handle Cells should be called
    bool p_run_cells;
    /// Flag if handle Cell Pairs should be called
    bool p_run_pairs;
    /// Reference to MoleculeContainer::soa
    SOA& p_soa;
};


#endif //KOMD_FORCEFUNCTOR_H
