//
// Created by alex on 7/31/24.
//

#include "MoleculeContainer.h"

#include "molecule/Molecule.h"
#include "Registry.h"

MoleculeContainer::MoleculeContainer() : m_data() {
    auto config = Registry::instance->configuration();

    // create cell structure
    math::d3 domain_size = (config->domainHigh - config->domainLow);
    math::ul3 num_cells_domain = math::ufloor(domain_size / config->cutoff) + 2;
    m_data.init(num_cells_domain);
    // set up cell bounds
    math::d3 low, high;
    low.z() = config->domainLow.z() - config->cutoff;
    for (uint64_t cz = 0; cz < num_cells_domain.z(); cz++) {
        high.z() = (cz == num_cells_domain.z() - 2) ? config->domainHigh.z() : (low.z() + config->cutoff);
        low.y() = config->domainLow.y() - config->cutoff;
        for (uint64_t cy = 0; cy < num_cells_domain.y(); cy++) {
            high.y() = (cy == num_cells_domain.y() - 2) ? config->domainHigh.y() : (low.y() + config->cutoff);
            low.x() = config->domainLow.x() - config->cutoff;
            for (uint64_t cx = 0; cx < num_cells_domain.x(); cx++) {
                high.x() = (cx == num_cells_domain.x() - 2) ? config->domainHigh.x() : (low.x() + config->cutoff);

                m_data(cx, cy, cz).setBounds(low, high);

                low.x() = high.x();
            }
            low.y() = high.y();
        }
        low.z() = high.z();
    }
}

void MoleculeContainer::addMolecule(const Molecule &molecule) {
    auto config = Registry::instance->configuration();

    // first find correct cell
    const math::d3 pos = molecule.getCenterOfMass();
    bool valid = false;
    math::ul3 cell_coord = findCell(pos, valid);
    if (!valid) throw std::runtime_error("Inserting Molecule at invalid position.");

    // insert into cell
    m_data(cell_coord).addMolecule(molecule);
}

void MoleculeContainer::updateContainer() {
    auto boundary = Registry::instance->boundary();

    // Move data back into AOS to move molecules
    writeSOA2AOS();

    // check all cells
    std::vector<Molecule> deletedMolecules;
    for (uint64_t z = 1; z < m_data.dims().z()-1; z++) {
        for (uint64_t y = 1; y < m_data.dims().y()-1; y++) {
            for (uint64_t x = 1; x < m_data.dims().x()-1; x++) {
                Cell& cell = m_data(x, y, z);
                auto& molecules = cell.molecules();
                for (auto it = molecules.begin(); it != molecules.end();) {
                    // check if molecule is within bounds of cell
                    const math::d3 mol_pos = it->getCenterOfMass();
                    if (cell.insideBounds(mol_pos)) {
                        ++it;
                        continue;
                    }

                    // not in cell -> must delete (and move optionally)
                    Molecule molecule = *it;
                    it = boundary->deleteMolecule(molecule);
                    //boundary->moveMolecule(molecule);
                    deletedMolecules.push_back(molecule);
                }
            }
        }
    }

    // update halo now, to not touch work twice
    boundary->updateHalo();

    // reinsert molecules scheduled for moving
    for (Molecule& molecule : deletedMolecules) boundary->moveMolecule(molecule);

    // reconstruct all broken SOAs
    constructSOAs();

    clearForces();
}

MoleculeContainer::Iterator MoleculeContainer::iterator(const IteratorType type, const IteratorRegion region) {
    math::ul3 low, high;
    if (region == DOMAIN) {
        low = {1, 1, 1};
        high = m_data.dims() - 2;
    }
    if (region == DOMAIN_HALO) {
        low = {0, 0, 0};
        high = m_data.dims() - 1;
    }

    if (type == MOLECULE) writeSOA2AOS();

    return {low, high, type == MOLECULE, m_data};
}

void MoleculeContainer::constructSOAs() {
    if (!Registry::instance->configuration()->enableSOA) return;

    for (uint64_t z = 0; z < m_data.dims().z(); z++) {
        for (uint64_t y = 0; y < m_data.dims().y(); y++) {
            for (uint64_t x = 0; x < m_data.dims().x(); x++) {
                m_data(x, y, z).constructSOA();
            }
        }
    }
}

void MoleculeContainer::writeSOA2AOS() {
    if (!Registry::instance->configuration()->enableSOA) return;

    for (uint64_t z = 0; z < m_data.dims().z(); z++) {
        for (uint64_t y = 0; y < m_data.dims().y(); y++) {
            for (uint64_t x = 0; x < m_data.dims().x(); x++) {
                m_data(x, y, z).writeSOA2AOS();
            }
        }
    }
}

math::ul3 MoleculeContainer::findCell(const math::d3 &pos, bool& valid) {
    auto config = Registry::instance->configuration();

    // first find correct cell
    math::ul3 cell_coord {0, 0, 0};
    valid = true;

    for (int dim = 0; dim < 3; dim++) {
        // check for invalid
        if (pos[dim] < config->domainLow[dim] - config->cutoff ||
            pos[dim] > config->domainHigh[dim] + config->cutoff) {
            valid = false;
            return cell_coord;
        }

        // get cell coord in this dimension
        if (pos[dim] < config->domainLow[dim]) cell_coord[dim] = 0;
        else if (pos[dim] > config->domainHigh[dim]) cell_coord[dim] = m_data.dims()[dim] - 1;
        else {
            uint64_t idx = static_cast<uint64_t>(std::floor((pos[dim] - config->domainLow[dim]) / config->cutoff)) + 1UL;
            if (idx >= m_data.dims()[dim] - 1) idx = m_data.dims()[dim] - 2;
            cell_coord[dim] = idx;
        }
    }

    return cell_coord;
}

Vec3D<Cell> &MoleculeContainer::getCells() {
    return m_data;
}

void MoleculeContainer::clearForces() {
    for (uint64_t z = 1; z < m_data.dims().z()-1; z++) {
        for (uint64_t y = 1; y < m_data.dims().y()-1; y++) {
            for (uint64_t x = 1; x < m_data.dims().x()-1; x++) {
                Cell& cell = m_data(x, y, z);
                cell.clearForces();
            }
        }
    }
}

MoleculeContainer::Iterator::Iterator(const math::ul3 &min, const math::ul3 &max, bool only_mol, Vec3D<Cell> &cells) :
m_cell_coord(min), m_cell_min(min), m_cell_max(max), m_only_molecule(only_mol),
m_site_idx(0), m_molecule_idx(0), m_visited_sites_cell(0), m_cells(cells), m_cell_has_soa(cells(min).validSOA()) {
    if (!cells(min).molecules().empty()) return;

    // cell is empty, find the first cell in iteration direction that is not empty
    findNextCell();
    m_cell_has_soa = cells(m_cell_coord).validSOA();
}

void MoleculeContainer::Iterator::operator++() {
    Cell& cell = m_cells(m_cell_coord);
    auto& cell_molecules = cell.molecules();

    // move to next site
    m_visited_sites_cell++;
    if (!m_only_molecule) {
        if (m_site_idx < cell_molecules[m_molecule_idx].getSites().size() - 1) {
            m_site_idx++;
            return;
        }
    }

    // move to next molecule
    if (m_molecule_idx < cell_molecules.size() - 1) {
        m_molecule_idx++;
        m_site_idx = 0;
        return;
    }

    // move to next cell
    findNextCell();
    m_cell_has_soa = m_cells(m_cell_coord).validSOA();
}

bool MoleculeContainer::Iterator::isValid() {
    // check if cell is invalid
    if (m_cell_coord.x() == m_cell_max.x() + 1 &&
        m_cell_coord.y() == m_cell_max.y() + 1 &&
        m_cell_coord.z() == m_cell_max.z() + 1) return false;

    Cell& cell = m_cells(m_cell_coord);
    auto& cell_molecules = cell.molecules();
    // cell is valid
    // check if current molecule has more sites
    if (!m_only_molecule) // only do this if we iterate over sites
        if (m_site_idx < cell_molecules[m_molecule_idx].getSites().size() - 1) return true;

    // no more sites
    // check if there is another molecule in this cell
    if (m_molecule_idx < cell_molecules.size() - 1) return true;

    // no more molecules in this cell
    // check if there is another cell with molecules
    bool init = false;
    for (uint64_t z = m_cell_coord.z(); z <= m_cell_max.z(); z++) {
        for (uint64_t y = init ? m_cell_min.y() : m_cell_coord.y(); y <= m_cell_max.y(); y++) {
            for (uint64_t x = init ? m_cell_min.x() : m_cell_coord.x(); x <= m_cell_max.x(); x++) {
                if (m_cells(x, y, z).molecules().empty()) continue;
                // found non empty cell
                return true;
            }
            init = true;
        }
    }

    // all checks failed
    return false;
}

math::d3 & MoleculeContainer::Iterator::f() {
    if (m_cell_has_soa) return m_cells(m_cell_coord).soa().f()[m_visited_sites_cell];
    return m_cells(m_cell_coord).molecules()[m_molecule_idx].getSites()[m_site_idx].f_arr();
}

math::d3 & MoleculeContainer::Iterator::r() {
    if (m_cell_has_soa) return m_cells(m_cell_coord).soa().r()[m_visited_sites_cell];
    return m_cells(m_cell_coord).molecules()[m_molecule_idx].getSites()[m_site_idx].r_arr();
}

math::d3 & MoleculeContainer::Iterator::v() {
    if (m_cell_has_soa) return m_cells(m_cell_coord).soa().v()[m_visited_sites_cell];
    return m_cells(m_cell_coord).molecules()[m_molecule_idx].getSites()[m_site_idx].v_arr();
}

double MoleculeContainer::Iterator::epsilon() {
    if (m_cell_has_soa) return m_cells(m_cell_coord).soa().epsilon()[m_visited_sites_cell];
    return m_cells(m_cell_coord).molecules()[m_molecule_idx].getSites()[m_site_idx].getEpsilon();
}

double MoleculeContainer::Iterator::sigma() {
    if (m_cell_has_soa) return m_cells(m_cell_coord).soa().sigma()[m_visited_sites_cell];
    return m_cells(m_cell_coord).molecules()[m_molecule_idx].getSites()[m_site_idx].getSigma();
}

double MoleculeContainer::Iterator::mass() {
    if (m_cell_has_soa) return m_cells(m_cell_coord).soa().mass()[m_visited_sites_cell];
    return m_cells(m_cell_coord).molecules()[m_molecule_idx].getSites()[m_site_idx].getMass();
}

uint64_t MoleculeContainer::Iterator::ID() {
    if (m_cell_has_soa) return m_cells(m_cell_coord).soa().id()[m_visited_sites_cell];
    return m_cells(m_cell_coord).molecules()[m_molecule_idx].ID();
}

Molecule & MoleculeContainer::Iterator::molecule() {
    return m_cells(m_cell_coord).molecules()[m_molecule_idx];
}

void MoleculeContainer::Iterator::findNextCell() {
    bool init = false;
    for (uint64_t z = m_cell_coord.z(); z <= m_cell_max.z(); z++) {
        for (uint64_t y = init ? m_cell_min.y() : m_cell_coord.y(); y <= m_cell_max.y(); y++) {
            for (uint64_t x = init ? m_cell_min.x() : (m_cell_coord.x()+1); x <= m_cell_max.x(); x++) {
                if (m_cells(x, y, z).molecules().empty()) continue;
                // found non empty cell
                m_cell_coord = {x, y, z};
                m_molecule_idx = 0;
                m_site_idx = 0;
                m_visited_sites_cell = 0;
                return;
            }
            init = true;
        }
    }
    // if loop is not successful cell_coord = max+1 => isValid becomes false
    m_cell_coord = m_cell_max + 1;
}
