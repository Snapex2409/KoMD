//
// Created by alex on 8/20/24.
//

#include <iostream>
#include "LinkedCells.h"
#include "molecule/Molecule.h"
#include "Registry.h"
#include "math/Geometry.h"

LinkedCells::LinkedCells() : MoleculeContainer(), m_data() {
    auto config = Registry::instance->configuration();
    const int num_halo_cells = 2;

    // create cell structure
    m_low = config->domainLow;
    m_high = config->domainHigh;
    m_dom_size = m_high - m_low;
    const math::ul3 num_cells_domain = math::ufloor(m_dom_size / config->cell_size);
    const math::ul3 num_cells = num_cells_domain + num_halo_cells;
    m_data.init(num_cells);
    m_cutoff = config->cell_size;

    // set up cell bounds
    math::d3 low, high;
    const double halo_size = m_cutoff;
    low.z() = m_low.z() - halo_size;
    for (uint64_t cz = 0; cz < num_cells.z(); cz++) {
        high.z() = (cz == num_cells.z() - 2) ? m_high.z() : (low.z() + m_cutoff);
        low.y() = m_low.y() - halo_size;
        for (uint64_t cy = 0; cy < num_cells.y(); cy++) {
            high.y() = (cy == num_cells.y() - 2) ? m_high.y() : (low.y() + m_cutoff);
            low.x() = m_low.x() - halo_size;
            for (uint64_t cx = 0; cx < num_cells.x(); cx++) {
                high.x() = (cx == num_cells.x() - 2) ? m_high.x() : (low.x() + m_cutoff);
                m_data(cx, cy, cz).setBounds(low, high);
                m_data(cx, cy, cz).setCoords({cx, cy, cz});

                low.x() = high.x();
            }
            low.y() = high.y();
        }
        low.z() = high.z();
    }
}

void LinkedCells::init() {
    uint64_t num_sites = 0;
    for (auto& molecule : p_molecules) num_sites += molecule.getSites().size();
    p_soa.createBuffers(num_sites);
    constructSOAs();
    for (auto c_it = iteratorCell(); c_it->isValid(); ++(*c_it)) c_it->cell().createIndexBuffers(num_sites);

    // write indices
    p_com = KW::vec_t<math::d3>("Center of Masses", p_soa.size());
    updateCOM();
    writeIndices();

    // init pair list buffers
    p_pair_list.init(num_sites);
}

void LinkedCells::updateContainer() {
    resetIndices();

    updateCOM();

    // apply periodic boundary kernel
    Kokkos::parallel_for("Periodic Bound", p_soa.size(), LinkedCells::Periodic_Kernel(p_com, p_soa.r(), m_low, m_high, m_dom_size));
    // apply force reset kernel
    Kokkos::parallel_for("Force Reset", p_soa.size(), LinkedCells::FReset_Kernel(p_soa.f()));
    Kokkos::fence("OneCell - update");

    updateCOM();
    writeIndices();

    updatePairList();
}

void LinkedCells::writeIndices() {
    const auto num_cells_domain = m_data.dims() - 2;

    uint64_t s_idx = 0;
    for (uint64_t m_idx = 0; m_idx < p_molecule_count; m_idx++) {
        Molecule& molecule = p_molecules[m_idx];
        // first find correct cell
        const math::d3 com_pos = p_com[s_idx];
        math::ul3 cell_coord;
        for (int dim = 0; dim < 3; dim++) {
            double fBin = (com_pos[dim] - m_low[dim]) / m_cutoff;
            cell_coord[dim] = static_cast<uint64_t>(std::clamp((int64_t) fBin, 0L, (int64_t) num_cells_domain[dim] - 1));
        }
        // offset cell coord by 1 as we have one layer of halo
        cell_coord += 1;

        for (auto& _ : molecule.getSites()) {
            m_data(cell_coord).addIndex(s_idx);
            s_idx++;
        }
    }

    // all domain cells are populated
    // now we generate halo indices
    for (auto c_it = iteratorCell(); c_it->isValid(); ++(*c_it)) {
        Cell& cell = c_it->cell();
        const auto check_low = cell.low() < m_low;
        const auto check_high = cell.high() > m_high;
        const auto is_halo = check_low - check_high; // will be 1 if is halo low, 0 if domain, -1 if is halo high

        if (is_halo == 0) continue;

        const auto src_coord = cell.coord() + is_halo * num_cells_domain;
        Cell& src_cell = m_data(src_coord);
        for (int idx = 0; idx < src_cell.getNumIndices(); idx++) {
            cell.addIndex(src_cell.indices()[idx]);
        }
    }
}

void LinkedCells::resetIndices() {
    for (auto c_it = iteratorCell(); c_it->isValid(); ++(*c_it)) {
        Cell& cell = c_it->cell();
        cell.resetIndices();
    }
}

void LinkedCells::updatePairList() {
    if (p_pair_list.requiresUpdate()) {
        p_pair_list.reset();


        // check how much memory will be required
        uint64_t num_pairs = 0;
        for (auto it = iteratorC08(); it->isValid(); ++(*it)) {
            Cell& cell0 = it->cell0();
            Cell& cell1 = it->cell1();
            num_pairs += cell0.getNumIndices() * cell1.getNumIndices();
        }

        for (auto it = iteratorCell(); it->isValid(); ++(*it)) {
            Cell& cell = it->cell();
            const auto check_low = cell.low() < m_low;
            const auto check_high = cell.high() > m_high;
            const auto is_halo = check_low - check_high; // will be 1 if is halo low, 0 if domain, -1 if is halo high
            if (is_halo != 0) continue;

            const auto n = cell.getNumIndices();
            num_pairs += n * (n-1) / 2;
        }
        p_pair_list.resize(num_pairs);


        //============================
        // populate pair list

        // create cell pair kernels
        std::vector<PairListPair_Kernel> pair_kernels;
        PairListPair_Kernel pair_kernel;
        pair_kernel.stored_cell_pairs = 0;

        for (auto it = iteratorC08(); it->isValid(); ++(*it)) {
            Cell& cell0 = it->cell0();
            Cell& cell1 = it->cell1();
            const math::d3 shift0 = it->getCell0Shift();
            const math::d3 shift1 = it->getCell1Shift();

            if (pair_kernel.stored_cell_pairs == PairListPair_Kernel::MAX_PAIRS) {
                pair_kernels.push_back(pair_kernel);
                pair_kernel = PairListPair_Kernel();
                pair_kernel.stored_cell_pairs = 0;
            }

            const int write_idx = pair_kernel.stored_cell_pairs;
            // cell data
            pair_kernel.cell_pairs[write_idx].indices0 = cell0.indices();
            pair_kernel.cell_pairs[write_idx].indices1 = cell1.indices();
            pair_kernel.cell_pairs[write_idx].shift0 = shift0;
            pair_kernel.cell_pairs[write_idx].shift1 = shift1;

            // meta data
            pair_kernel.pair_counts[write_idx] = cell0.getNumIndices() * cell1.getNumIndices();
            int count_acc = 0;
            if (write_idx > 0) count_acc = pair_kernel.pair_counts_accumulated[write_idx-1] + pair_kernel.pair_counts[write_idx-1];
            pair_kernel.pair_counts_accumulated[write_idx] = count_acc;
            pair_kernel.pair_dims[write_idx][0] = cell0.getNumIndices();
            pair_kernel.pair_dims[write_idx][1] = cell1.getNumIndices();
            pair_kernel.stored_cell_pairs++;
        }
        if (pair_kernel.stored_cell_pairs != 0) pair_kernels.push_back(pair_kernel);

        // create single cell kernels
        std::vector<PairListSingle_Kernel> single_kernels;
        PairListSingle_Kernel single_kernel;
        single_kernel.stored_cells = 0;

        for (auto it = iteratorCell(); it->isValid(); ++(*it)) {
            Cell& cell = it->cell();
            const auto check_low = cell.low() < m_low;
            const auto check_high = cell.high() > m_high;
            const auto is_halo = check_low - check_high; // will be 1 if is halo low, 0 if domain, -1 if is halo high
            if (is_halo != 0) continue;

            if (single_kernel.stored_cells == PairListSingle_Kernel::MAX_PAIRS) {
                single_kernels.push_back(single_kernel);
                single_kernel = PairListSingle_Kernel();
                single_kernel.stored_cells = 0;
            }

            const int write_idx = single_kernel.stored_cells;
            // cell data
            single_kernel.cells[write_idx].indices = cell.indices();

            // meta data
            const int n = cell.getNumIndices();
            single_kernel.counts[write_idx] = n * (n-1) / 2;
            int count_acc = 0;
            if (write_idx > 0) count_acc = single_kernel.counts_accumulated[write_idx-1] + single_kernel.counts[write_idx-1];
            single_kernel.counts_accumulated[write_idx] = count_acc;
            single_kernel.stored_cells++;
        }
        if (single_kernel.stored_cells != 0) single_kernels.push_back(single_kernel);

        // set up globals and offsets
        KW::vec_t<uint64_t> global_pair_idx = KW::vec_t<uint64_t>("PLU global idx", 1);
        global_pair_idx(0) = 0;
        uint64_t write_offset = 0;
        PairList_Globals globals { p_pair_list.getPairs(), p_pair_list.getOffsets(), global_pair_idx };
        for (auto& kernel : pair_kernels) {
            kernel.globals = globals;
            kernel.write_offset = write_offset;
            write_offset += kernel.pair_counts_accumulated[kernel.stored_cell_pairs-1] + kernel.pair_counts[kernel.stored_cell_pairs-1];
        }
        for (auto& kernel : single_kernels) {
            kernel.globals = globals;
            kernel.write_offset = write_offset;
            write_offset += kernel.counts_accumulated[kernel.stored_cells-1] + kernel.counts[kernel.stored_cells-1];
        }

        // run kernels
        for (auto& kernel : pair_kernels) Kokkos::parallel_for("PLU - Pair", kernel.pair_counts_accumulated[kernel.stored_cell_pairs-1] + kernel.pair_counts[kernel.stored_cell_pairs-1], kernel);
        for (auto& kernel : single_kernels) Kokkos::parallel_for("PLU - Single", kernel.counts_accumulated[kernel.stored_cells-1] + kernel.counts[kernel.stored_cells-1], kernel);
        Kokkos::fence("PLU - End");
    }

    p_pair_list.step();
}

void LinkedCells::FReset_Kernel::operator()(int idx) const {
    f[idx] = math::d3 {0, 0, 0};
}

void LinkedCells::Periodic_Kernel::operator()(int idx) const {
    const math::d3 pos = com[idx];

    math::d3 offset {0, 0, 0};
    for (int dim = 0; dim < 3; dim++) {
        if (pos[dim] < low[dim]) offset[dim] = domain_size[dim];
        if (pos[dim] > high[dim]) offset[dim] = -domain_size[dim];
    }

    r[idx] += offset;
}

//===========================================================================
// Iterator
//===========================================================================

std::unique_ptr<MoleculeContainer::CellIterator> LinkedCells::iteratorCell() {
    math::ul3 low, high;
    low = {0, 0, 0};
    high = m_data.dims()-1;

    return std::make_unique<LCCellIterator>(low, high, m_data);
}

std::unique_ptr<MoleculeContainer::CellPairIterator> LinkedCells::iteratorC08() {
    return std::make_unique<LCC08Iterator>(math::ul3{0,0,0}, m_data.dims() - 1, m_data, m_dom_size);
}


LinkedCells::LCCellIterator::LCCellIterator(const math::ul3 &min, const math::ul3 &max, Vec3D<Cell> &cells) :
        m_cell_coord(min), m_cell_min(min), m_cell_max(max), m_cells(cells) { }

void LinkedCells::LCCellIterator::operator++() {
    if (m_cell_coord.x() < m_cell_max.x()) { m_cell_coord.x() += 1; return; }
    m_cell_coord.x() = m_cell_min.x();
    if (m_cell_coord.y() < m_cell_max.y()) { m_cell_coord.y() += 1; return; }
    m_cell_coord.y() = m_cell_min.y();
    if (m_cell_coord.z() < m_cell_max.z()) { m_cell_coord.z() += 1; return; }
    m_cell_coord = m_cell_max + 1;
}

bool LinkedCells::LCCellIterator::isValid() const {
    return m_cell_coord != (m_cell_max + 1);
}

Cell &LinkedCells::LCCellIterator::cell() {
    return m_cells(m_cell_coord);
}

LinkedCells::LCC08Iterator::LCC08Iterator(const math::ul3 &min, const math::ul3 &max, Vec3D<Cell> &cells, const math::d3& dom_size) :
        m_cell_coord(min), m_cell_min(min), m_cell_max(max), m_cells(cells), m_offset_idx(0), m_color(0), m_color_switched(false), m_dom_size(dom_size), m_cell_dims(cells.dims()) {
    checkState();
}

void LinkedCells::LCC08Iterator::operator++() {
    m_color_switched = false;
    if (m_offset_idx < pair_offsets.size() - 1) {
        m_offset_idx++;
        return checkState();
    }
    m_offset_idx = 0;

    if (m_cell_coord.x() < m_cell_max.x()-2) { m_cell_coord.x() += 2; return checkState(); }
    m_cell_coord.x() = m_cell_min.x();
    if (m_cell_coord.y() < m_cell_max.y()-2) { m_cell_coord.y() += 2; return checkState(); }
    m_cell_coord.y() = m_cell_min.y();
    if (m_cell_coord.z() < m_cell_max.z()-2) { m_cell_coord.z() += 2; return checkState(); }

    m_color++;
    const math::i3 begin {m_color & 0b1, (m_color & 0b10) >> 1, (m_color & 0b100) >> 2};
    m_cell_min = begin;
    m_cell_coord = begin;
    m_color_switched = true;
    return checkState();
}


void LinkedCells::LCC08Iterator::checkState() {
    if (!isValid()) return;
    if (getCell0Shift() != 0 && getCell1Shift() != 0) return this->operator++();
}

bool LinkedCells::LCC08Iterator::isValid() const {
    return m_color < 8;
}

Cell &LinkedCells::LCC08Iterator::cell0() {
    const math::ul3 coord = m_cell_coord + pair_offsets[m_offset_idx].first;
    return m_cells(coord);
}

Cell &LinkedCells::LCC08Iterator::cell1() {
    const math::ul3 coord = m_cell_coord + pair_offsets[m_offset_idx].second;
    return m_cells(coord);
}

bool LinkedCells::LCC08Iterator::colorSwitched() {
    return m_color_switched;
}

math::d3 LinkedCells::LCC08Iterator::getCell0Shift() {
    const math::ul3 coord = m_cell_coord + pair_offsets[m_offset_idx].first;
    return m_dom_size * (coord >= m_cell_max) - m_dom_size * (coord <= math::ul3{0, 0, 0});
}

math::d3 LinkedCells::LCC08Iterator::getCell1Shift() {
    const math::ul3 coord = m_cell_coord + pair_offsets[m_offset_idx].second;
    return m_dom_size * (coord >= m_cell_max) - m_dom_size * (coord <= math::ul3{0, 0, 0});
}
