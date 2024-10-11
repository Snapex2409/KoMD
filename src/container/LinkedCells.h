//
// Created by alex on 8/20/24.
//

#ifndef KOMD_LINKEDCELLS_H
#define KOMD_LINKEDCELLS_H

#include <vector>
#include <memory>
#include <set>

#include "Vec3D.h"
#include "Cell.h"
#include "molecule/Molecule.h"
#include "MoleculeContainer.h"

class LinkedCells : public MoleculeContainer {
public:
    LinkedCells();
    void updateContainer() override;
    void init() override;

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

    /**
     * Does not guarantee overlapping writes between different cells, but reduces amount due to halo implementation
     * */
    class LCC08Iterator : public CellPairIterator {
    public:
        LCC08Iterator(const math::ul3& min, const math::ul3& max, Vec3D<Cell>& cells, const math::d3& dom_size);
        ~LCC08Iterator() override = default;
        void operator++() override;
        bool isValid() const override;
        Cell& cell0() override;
        Cell& cell1() override;
        bool colorSwitched() override;
        math::d3 getCell0Shift() override;
        math::d3 getCell1Shift() override;
    private:
        void checkState();

        math::ul3 m_cell_coord;
        math::ul3 m_cell_min;
        const math::ul3 m_cell_max;
        Vec3D<Cell>& m_cells;
        int m_offset_idx;
        int m_color;
        bool m_color_switched;
        const math::d3 m_dom_size;
        const math::ul3 m_cell_dims;

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
     * Only to be used on CPU side
     * */
    std::unique_ptr<CellIterator> iteratorCell() override;
    /**
     * Only to be used on CPU side
     * */
    std::unique_ptr<CellPairIterator> iteratorC08() override;

    struct FReset_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> f;
    };

    struct Periodic_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> com;
        KW::vec_t<math::d3> r;
        const math::d3 low;
        const math::d3 high;
        const math::d3 domain_size;
    };

    struct PairListPair_Kernel {
        KOKKOS_FUNCTION void operator()(int idx0, int idx1) const {
            uint64_t updated_pos = Kokkos::atomic_add_fetch(&num_pairs(0), 1);
            const uint64_t local_pos = updated_pos - 1;

            pairs(local_pos, 0) = indices0(idx0);
            pairs(local_pos, 1) = indices1(idx1);
            offsets(local_pos, 0) = shift0;
            offsets(local_pos, 1) = shift1;
        }
        const math::d3 shift0;
        const math::d3 shift1;
        KW::vec_t<uint64_t> indices0;
        KW::vec_t<uint64_t> indices1;
        KW::nvec_t<int, 2> pairs;
        KW::nvec_t<math::d3, 2> offsets;
        KW::vec_t<uint64_t> num_pairs;
    };

    struct PairListSingle_Kernel {
        KOKKOS_FUNCTION void operator()(int idx0, int idx1) const {
            if (idx1 <= idx0) return;
            uint64_t updated_pos = Kokkos::atomic_add_fetch(&num_pairs(0), 1);
            const uint64_t local_pos = updated_pos - 1;

            pairs(local_pos, 0) = indices(idx0);
            pairs(local_pos, 1) = indices(idx1);
            offsets(local_pos, 0) = {0, 0, 0};
            offsets(local_pos, 1) = {0, 0, 0};
        }
        KW::vec_t<uint64_t> indices;
        KW::nvec_t<int, 2> pairs;
        KW::nvec_t<math::d3, 2> offsets;
        KW::vec_t<uint64_t> num_pairs;
    };
private:
    /**
     * Writes all soa indices into cell buffers
     * */
    void writeIndices();

    /**
     * Resets all cell index buffers (does not reallocate memory, only sets counters to 0)
     * */
    void resetIndices();

    /**
     * Updates Pair List based on newly created Indices
     * */
    void updatePairList();

    /// 3d cell buffer
    Vec3D<Cell> m_data;
    /// domain low
    math::d3 m_low;
    /// domain high
    math::d3 m_high;
    /// domain size
    math::d3 m_dom_size;
    /// cutoff
    double m_cutoff;
};


#endif //KOMD_LINKEDCELLS_H
