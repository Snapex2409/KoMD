//
// Created by alex on 7/31/24.
//

#include "ForceFunctor.h"
#include "Registry.h"

ForceFunctor::ForceFunctor() : p_use_soa(Registry::instance->configuration()->enableSOA), p_run_cells(true), p_run_pairs(true) {}

void ForceFunctor::operator()() {
    if (p_run_cells) iterateCells();
    if (p_run_pairs) iterateCellPairs();
}

void ForceFunctor::iterateCells() {
    auto& cells = Registry::instance->moleculeContainer()->getCells();
    const math::ul3 cell_dims = cells.dims();

    for (uint64_t z = 1; z < cell_dims.z()-1; z++) {
        for (uint64_t y = 1; y < cell_dims.y()-1; y++) {
            for (uint64_t x = 1; x < cell_dims.x()-1; x++) {
                Cell& cell = cells[x, y, z];
                handleCell(cell);
            }
        }
    }
}

void ForceFunctor::iterateCellPairs() {
    auto& cells = Registry::instance->moleculeContainer()->getCells();
    const math::ul3 cell_dims = cells.dims();

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

    //we are implementing C08 traversal here
    for (uint64_t col = 0; col < 8; col++) {
        const math::ul3 begin {col & 0b1, (col & 0b10) >> 1, (col & 0b100) >> 2};
        const math::ul3 end = cell_dims - 1;

        for (uint64_t z = begin.z(); z < end.z(); z += 2) {
            for (uint64_t y = begin.y(); y < end.y(); y += 2) {
                for (uint64_t x = begin.x(); x < end.x(); x += 2) {
                    const math::ul3 base_coord {x, y, z};
                    for (auto& [offset_first, offset_second] : pair_offsets) {
                        const math::ul3 coordFirst = base_coord + offset_first;
                        const math::ul3 coordSecond = base_coord + offset_second;

                        Cell& cellFirst = cells[coordFirst];
                        Cell& cellSecond = cells[coordSecond];
                        handleCellPair(cellFirst, cellSecond);
                    }
                }
            }
        }
    }
}
