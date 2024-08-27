//
// Created by alex on 8/20/24.
//

#include "LinkedCells.h"
#include "molecule/Molecule.h"
#include "Registry.h"

LinkedCells::LinkedCells() : MoleculeContainer(), m_data() {
    auto config = Registry::instance->configuration();

    // create cell structure
    m_low = config->domainLow;
    m_high = config->domainHigh;
    m_dom_size = m_high - m_low;
    math::ul3 num_cells_domain = math::ufloor(m_dom_size / config->cutoff);
    m_data.init(num_cells_domain);
    m_cutoff = config->cutoff;

    // set up cell bounds
    math::d3 low, high;
    low.z() = m_low.z();
    for (uint64_t cz = 0; cz < num_cells_domain.z(); cz++) {
        high.z() = (cz == num_cells_domain.z() - 1) ? m_high.z() : (low.z() + m_cutoff);
        low.y() = m_low.y();
        for (uint64_t cy = 0; cy < num_cells_domain.y(); cy++) {
            high.y() = (cy == num_cells_domain.y() - 1) ? m_high.y() : (low.y() + m_cutoff);
            low.x() = m_low.x();
            for (uint64_t cx = 0; cx < num_cells_domain.x(); cx++) {
                high.x() = (cx == num_cells_domain.x() - 1) ? m_high.x() : (low.x() + m_cutoff);
                m_data(cx, cy, cz).setBounds(low, high);
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
    updateCOM();
    writeIndices();
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
}

void LinkedCells::writeIndices() {
    uint64_t s_idx = 0;
    for (uint64_t m_idx = 0; m_idx < p_molecule_count; m_idx++) {
        Molecule& molecule = p_molecules[m_idx];
        // first find correct cell
        const math::d3 com_pos = p_com[s_idx];
        math::ul3 cell_coord;
        for (int dim = 0; dim < 3; dim++) {
            double fBin = (com_pos[dim] - m_low[dim]) / m_cutoff;
            cell_coord[dim] = static_cast<uint64_t>(std::clamp((int64_t) fBin, 0L, (int64_t) m_data.dims()[dim]-1));
        }

        for (auto& _ : molecule.getSites()) {
            m_data(cell_coord).addIndex(s_idx);
            s_idx++;
        }
    }
}

void LinkedCells::resetIndices() {
    for (auto c_it = iteratorCell(); c_it->isValid(); ++(*c_it)) {
        Cell& cell = c_it->cell();
        cell.resetIndices();
    }
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
        m_cell_coord(min), m_cell_min(min), m_cell_max(max), m_cells(cells), m_offset_idx(0), m_color(0), m_color_switched(false), m_dom_size(dom_size), m_cell_dims(cells.dims()) { }

void LinkedCells::LCC08Iterator::operator++() {
    m_color_switched = false;
    if (m_offset_idx < pair_offsets.size() - 1) { m_offset_idx++; return; }
    m_offset_idx = 0;

    // -1 is intended -> allow for periodic bound cell shift, is handled during cell access
    if (m_cell_coord.x() < m_cell_max.x()-1) { m_cell_coord.x() += 2; return; }
    m_cell_coord.x() = m_cell_min.x();
    if (m_cell_coord.y() < m_cell_max.y()-1) { m_cell_coord.y() += 2; return; }
    m_cell_coord.y() = m_cell_min.y();
    if (m_cell_coord.z() < m_cell_max.z()-1) { m_cell_coord.z() += 2; return; }

    m_color++;
    const math::i3 begin {m_color & 0b1, (m_color & 0b10) >> 1, (m_color & 0b100) >> 2};
    m_cell_min = begin;
    m_cell_coord = begin;
    m_color_switched = true;
}

bool LinkedCells::LCC08Iterator::isValid() const {
    return m_color < 8;
}

Cell &LinkedCells::LCC08Iterator::cell0() {
    const math::ul3 coord = m_cell_coord + pair_offsets[m_offset_idx].first;
    return m_cells(coord);
}

Cell &LinkedCells::LCC08Iterator::cell1() {
    math::ul3 coord = m_cell_coord + pair_offsets[m_offset_idx].second;
    const math::ul3 required_shift_dims = coord > m_cell_max;
    if (required_shift_dims != 0) coord -= required_shift_dims * m_cell_dims;
    return m_cells(coord);
}

bool LinkedCells::LCC08Iterator::colorSwitched() {
    return m_color_switched;
}

math::d3 LinkedCells::LCC08Iterator::getCell1Shift() {
    const math::ul3 coord = m_cell_coord + pair_offsets[m_offset_idx].second;
    const math::ul3 required_shift_dims = coord > m_cell_max;
    return m_dom_size * required_shift_dims;
}
