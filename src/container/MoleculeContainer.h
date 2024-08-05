//
// Created by alex on 7/31/24.
//

#ifndef KOMD_MOLECULECONTAINER_H
#define KOMD_MOLECULECONTAINER_H

#include <vector>
#include <memory>

#include "Vec3D.h"
#include "Cell.h"
#include "molecule/Molecule.h"

class MoleculeContainer {
public:
    MoleculeContainer();
    void addMolecule(const Molecule& molecule);
    void updateContainer();
    Vec3D<Cell>& getCells();
    /**
     * Returns a valid cell coordinate if pos is within simulation bounds, else returned result is invalid and valid is set to false
     * */
    math::ul3 findCell(const math::d3& pos, bool& valid);
private:
    Vec3D<Cell> m_data;
    void constructSOAs();
    void writeSOA2AOS();
};


#endif //KOMD_MOLECULECONTAINER_H
