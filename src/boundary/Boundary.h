//
// Created by alex on 8/1/24.
//

#ifndef KOMD_BOUNDARY_H
#define KOMD_BOUNDARY_H

#include <functional>
#include "math/Array.h"

class Cell;
class Molecule;

class Boundary {
public:
    /**
     * Called once before simulation start, to set up Halo correctly.\n
     * Must be called after Molecule container is constructed and all desired molecules were inserted.
     * */
    void setup();

    /**
     * Moves a molecule into the correct cell. If no such cell exists, molecule will not be added.\n
     * Will check if domain bounds are reached. If so, then wrapping will be handled by this method.
     * */
    void moveMolecule(Molecule& molecule);

    /**
     * Moves all Halo molecules to the correct position
     * */
    void updateHalo();

    /**
     * Deletes this molecule and all copies from all cells
     * */
    std::vector<Molecule>::iterator deleteMolecule(Molecule& molecule);

    /**
     * Creates HaloMolecules for the given molecule, if it is in a boundary cell.
     * Otherwise, nothing happens.
     * @param cell_coord coordinate of cell, in which molecule is
     * @param domain_size size of active domain
     * */
    void createHaloMolecules(Molecule& molecule, const math::ul3& cell_coord, const math::d3& domain_size);

    /**
     * Patches HaloMolecules for the given molecule, if it is in a boundary cell. (Fixes the position)
     * Otherwise, nothing happens.
     * @param cell_coord coordinate of cell, in which molecule is
     * @param domain_size size of active domain
     * */
    void updateHaloMolecules(Molecule& molecule, const math::ul3& cell_coord, const math::d3& domain_size);

    /**
     * Loops over all boundary cells while executing fun. Provides the current cell and it coordinate as parameters.
     * */
    void loopOverBoundary(std::function<void(Cell&, const math::ul3&)> fun);
private:

};


#endif //KOMD_BOUNDARY_H
