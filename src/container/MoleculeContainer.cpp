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

                m_data[cx, cy, cz].setBounds(low, high);

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
    m_data[cell_coord].addMolecule(molecule);
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
                Cell& cell = m_data[x, y, z];
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

void MoleculeContainer::constructSOAs() {
    if (!Registry::instance->configuration()->enableSOA) return;

    for (uint64_t z = 0; z < m_data.dims().z(); z++) {
        for (uint64_t y = 0; y < m_data.dims().y(); y++) {
            for (uint64_t x = 0; x < m_data.dims().x(); x++) {
                m_data[x, y, z].constructSOA();
            }
        }
    }
}

void MoleculeContainer::writeSOA2AOS() {
    if (!Registry::instance->configuration()->enableSOA) return;

    for (uint64_t z = 0; z < m_data.dims().z(); z++) {
        for (uint64_t y = 0; y < m_data.dims().y(); y++) {
            for (uint64_t x = 0; x < m_data.dims().x(); x++) {
                m_data[x, y, z].writeSOA2AOS();
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
                Cell& cell = m_data[x, y, z];
                cell.clearForces();
            }
        }
    }
}
