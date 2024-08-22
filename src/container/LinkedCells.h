//
// Created by alex on 8/20/24.
//

#ifndef KOMD_LINKEDCELLS_H
#define KOMD_LINKEDCELLS_H

#include <vector>
#include <memory>

#include "Vec3D.h"
#include "Cell.h"
#include "molecule/Molecule.h"
#include "MoleculeContainer.h"
#include "boundary/LCBoundary.h"

class LinkedCells : public MoleculeContainer {
public:
    LinkedCells();
    void addMolecule(const Molecule& molecule) override;
    void updateContainer() override;
    void init() override;
    Vec3D<Cell>& getCells();
    math::ul3 findCell(const math::d3& pos, bool& valid);
    void getCenterOfMassPositions(SOA::vec_t<math::d3>& buffer) override;
    void writeSOA2AOS() override;

    class LCIterator : public Iterator {
    public:
        LCIterator(const math::ul3& min, const math::ul3& max, bool only_mol, Vec3D<Cell>& cells);
        ~LCIterator() override = default;
        void operator++() override;
        bool isValid() const override;

        math::d3& f() override;
        math::d3& r() override;
        math::d3& v() override;
        double epsilon() override;
        double sigma() override;
        double mass() override;
        uint64_t ID() override;

        Molecule& molecule() override;
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

    class LCCellIterator : public CellIterator {
    public:
        LCCellIterator(const math::ul3& min, const math::ul3& max, Vec3D<Cell>& cells);
        ~LCCellIterator() override = default;
        void operator++() override;
        bool isValid() const override;
        Cell& cell() override;
    private:
        math::ul3 m_cell_coord;
        const math::ul3 m_cell_min;
        const math::ul3 m_cell_max;
        Vec3D<Cell>& m_cells;
    };

    class LCC08Iterator : public CellPairIterator {
    public:
        LCC08Iterator(const math::ul3& min, const math::ul3& max, Vec3D<Cell>& cells);
        ~LCC08Iterator() override = default;
        void operator++() override;
        bool isValid() const override;
        Cell& cell0() override;
        Cell& cell1() override;
        bool colorSwitched() override;
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
    std::unique_ptr<Iterator> iterator(const IteratorType type, const IteratorRegion region) override;
    /**
     * Only to be used on CPU side
     * */
    std::unique_ptr<CellIterator> iteratorCell(const IteratorRegion region) override;
    /**
     * Only to be used on CPU side
     * */
    std::unique_ptr<CellPairIterator> iteratorC08() override;
private:
    /// 3d cell buffer
    Vec3D<Cell> m_data;
    void constructSOAs() override;
    void constructSOABuffers();
    void clearForces() override;
    LCBoundary m_boundary;
};


#endif //KOMD_LINKEDCELLS_H
