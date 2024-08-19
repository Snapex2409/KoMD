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
    uint64_t getNumMolecules() { return m_molecule_count; }
    /**
     * Returns a valid cell coordinate if pos is within simulation bounds, else returned result is invalid and valid is set to false
     * */
    math::ul3 findCell(const math::d3& pos, bool& valid);
    void gatherMolecules(std::vector<std::reference_wrapper<Molecule>>& buffer);
    void getCenterOfMassPositions(Kokkos::View<math::d3*, Kokkos::SharedSpace>& buffer);
    void writeSOA2AOS();

    enum IteratorType {
        MOLECULE,
        SITE
    };

    enum IteratorRegion {
        DOMAIN,
        DOMAIN_HALO
    };

    /**
     * Only to be used for IO functionality on CPU side
     * */
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
    };

    struct CellIterator {
    public:
        CellIterator(const math::ul3& min, const math::ul3& max, Vec3D<Cell>& cells);
        void operator++();
        bool isValid();
        Cell& cell();
    private:
        math::ul3 m_cell_coord;
        const math::ul3 m_cell_min;
        const math::ul3 m_cell_max;
        Vec3D<Cell>& m_cells;
    };

    struct C08Iterator {
    public:
        C08Iterator(const math::ul3& min, const math::ul3& max, Vec3D<Cell>& cells);
        void operator++();
        bool isValid() const;
        Cell& cell0();
        Cell& cell1();
        bool colorSwitched();
    private:
        math::ul3 m_cell_coord;
        math::ul3 m_cell_min;
        const math::ul3 m_cell_max;
        Vec3D<Cell>& m_cells;
        int m_offset_idx;
        int m_color;
        bool m_color_switched;

        static constexpr math::ul3 c_o   {0, 0, 0};
        static constexpr math::ul3 c_x   {1, 0, 0};
        static constexpr math::ul3 c_y   {0, 1, 0};
        static constexpr math::ul3 c_z   {0, 0, 1};
        static constexpr math::ul3 c_xy  {1, 1, 0};
        static constexpr math::ul3 c_yz  {0, 1, 1};
        static constexpr math::ul3 c_xz  {1, 0, 1};
        static constexpr math::ul3 c_xyz {1, 1, 1};

        using pair_t = std::pair<math::ul3, math::ul3>;
        static constexpr std::array<pair_t, 13> pair_offsets {pair_t {c_o, c_y},
                                                              pair_t {c_y, c_z},
                                                              pair_t {c_o, c_z},
                                                              pair_t {c_o, c_yz},
                                                              pair_t {c_x, c_yz},
                                                              pair_t {c_x, c_y},
                                                              pair_t {c_x, c_z},
                                                              pair_t {c_o, c_x},
                                                              pair_t {c_o, c_xy},
                                                              pair_t {c_xy, c_z},
                                                              pair_t {c_y, c_xz},
                                                              pair_t {c_o, c_xz},
                                                              pair_t {c_o, c_xyz}};
    };

    /**
     * Only to be used for IO functionality on CPU side
     * */
    Iterator iterator(const IteratorType type, const IteratorRegion region);
    CellIterator iteratorCell(const IteratorRegion region);
    C08Iterator iteratorC08();
private:
    Vec3D<Cell> m_data;
    void constructSOAs();
    void clearForces();
    uint64_t m_molecule_count;
};


#endif //KOMD_MOLECULECONTAINER_H
