//
// Created by alex on 7/31/24.
//

#ifndef KOMD_FORCEFUNCTOR_H
#define KOMD_FORCEFUNCTOR_H

class Cell;

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
    virtual void handleCellPair(Cell& cell0, Cell& cell1) = 0;

    /// Flag if soa should be used
    bool p_use_soa;
    /// Flag if handle Cells should be called
    bool p_run_cells;
    /// Flag if handle Cell Pairs should be called
    bool p_run_pairs;
private:
    void iterateCells();
    void iterateCellPairs();
};


#endif //KOMD_FORCEFUNCTOR_H
