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
    void writeSOA2AOS();

    struct Iterator {
    public:
        Iterator(const math::ul3& min, const math::ul3& max, bool only_mol, Vec3D<Cell>& cells);
        void operator++();
        bool isValid();

        math::d3& f();
        math::d3& r();
        math::d3& v();
        double epsilon();
        double sigma();
        double mass();
        uint64_t ID();

        Molecule& molecule();
    private:
        void findNextCell();

        math::ul3 m_cell_coord;
        const math::ul3 m_cell_min;
        const math::ul3 m_cell_max;
        const bool m_only_molecule;
        uint64_t m_site_idx;
        uint64_t m_molecule_idx;
        uint64_t m_visited_sites_cell;
        Vec3D<Cell>& m_cells;
        bool m_cell_has_soa;
    };

    enum IteratorType {
        MOLECULE,
        SITE
    };

    enum IteratorRegion {
        DOMAIN,
        DOMAIN_HALO
    };

    Iterator iterator(const IteratorType type, const IteratorRegion region);
private:
    Vec3D<Cell> m_data;
    void constructSOAs();
    void clearForces();
};


#endif //KOMD_MOLECULECONTAINER_H
