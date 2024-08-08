//
// Created by alex on 8/1/24.
//

#include "Boundary.h"

#include "molecule/Molecule.h"
#include "container/Cell.h"
#include "Registry.h"

void Boundary::setup() {
    // We want to find all molecules close to the boundary region and make copies to all sides
    auto config = Registry::instance->configuration();
    const math::d3 domain_size = config->domainHigh - config->domainLow;
    loopOverBoundary([&](Cell& cell, const math::ul3& cell_coord) -> void {
        for (uint64_t m_idx = 0; m_idx < cell.molecules().size(); m_idx++) {
            Molecule& molecule = cell.molecules()[m_idx];
            createHaloMolecules(molecule, cell, m_idx, cell_coord, domain_size);
        }
    });
}

void Boundary::updateHalo() {
    auto config = Registry::instance->configuration();
    const math::d3 domain_size = config->domainHigh - config->domainLow;
    loopOverBoundary([&](Cell& cell, const math::ul3& cell_coord) -> void {
        for (Molecule& molecule : cell.molecules()) {
            updateHaloMolecules(molecule, cell_coord, domain_size);
        }
    });
}

void Boundary::moveMolecule(Molecule& molecule) {
    auto config = Registry::instance->configuration();
    auto container = Registry::instance->moleculeContainer();

    // check if molecule is still within simulation bounds
    const math::d3 mol_pos = molecule.getCenterOfMass();
    for (int dim = 0; dim < 3; dim++) {
        if (mol_pos[dim] < config->domainLow[dim] - config->cutoff ||
            mol_pos[dim] > config->domainHigh[dim] + config->cutoff) {
            return;
        }
    }

    // handle wrapping, if necessary
    math::d3 offset {0, 0, 0};
    const math::d3 domain_size = config->domainHigh - config->domainLow;
    for (int dim = 0; dim < 3; dim++) {
        if (mol_pos[dim] < config->domainLow[dim]) offset[dim] = domain_size[dim];
        if (mol_pos[dim] > config->domainHigh[dim]) offset[dim] = -domain_size[dim];
    }
    molecule.moveBy(offset);

    // insert molecule based on new position
    bool valid = false;
    const math::d3 new_pos = molecule.getCenterOfMass();
    const math::ul3 cell_coord = container->findCell(new_pos, valid);
    Cell& cell = container->getCells()[cell_coord];
    cell.addMolecule(molecule);
    cell.invalidateSOA();

    // add halo molecules if needed
    createHaloMolecules(cell.molecules().back(), cell, cell.molecules().size()-1, cell_coord, domain_size);
}

std::vector<Molecule>::iterator Boundary::deleteMolecule(Molecule &molecule) {
    // remove all linked molecules
    for (auto& [cell_ref, _] : molecule.getCopies()) {
        Cell& cell = cell_ref.get();
        Molecule& mol = *cell.findBy(molecule.ID());
        cell.removeMolecule(mol.ID());
        cell.invalidateSOA();
    }
    molecule.getCopies().clear();

    // remove main molecule
    molecule.getCell().invalidateSOA();
    return molecule.getCell().removeMolecule(molecule.ID());
}

void Boundary::createHaloMolecules(Molecule& molecule, Cell& cell, uint64_t idx, const math::ul3& cell_coord, const math::d3& domain_size) {
    auto container = Registry::instance->moleculeContainer();

    const math::ul3 cell_dims = container->getCells().dims();
    math::d3 is_bound {0, 0, 0};
    if (cell_coord.x() == 1) is_bound.x() = 1;
    if (cell_coord.x() == cell_dims.x() - 2) is_bound.x() = -1;
    if (cell_coord.y() == 1) is_bound.y() = 1;
    if (cell_coord.y() == cell_dims.y() - 2) is_bound.y() = -1;
    if (cell_coord.z() == 1) is_bound.z() = 1;
    if (cell_coord.z() == cell_dims.z() - 2) is_bound.z() = -1;

    for (int shiftZ = 0; shiftZ <= ((is_bound.z() != 0) ? 1 : 0); shiftZ++) {
        for (int shiftY = 0; shiftY <= ((is_bound.y() != 0) ? 1 : 0); shiftY++) {
            for (int shiftX = 0; shiftX <= ((is_bound.x() != 0) ? 1 : 0); shiftX++) {
                if (shiftX == 0 && shiftY == 0 && shiftZ == 0) continue;
                const math::i3 shiftVec {shiftX, shiftY, shiftZ};
                const math::d3 halo_offset = domain_size * is_bound * shiftVec;

                Molecule halo_copy = molecule;
                halo_copy.moveBy(halo_offset);
                halo_copy.setParent(cell);

                const math::ul3 halo_cell_coord = cell_coord + (is_bound * shiftVec) * (cell_dims - 2);
                Cell& halo_cell = container->getCells()[halo_cell_coord];
                halo_cell.addMolecule(halo_copy);
                molecule.registerCopy(halo_cell, shiftVec);
                halo_cell.invalidateSOA();
            }
        }
    }
}

void Boundary::updateHaloMolecules(Molecule &molecule, const math::ul3 &cell_coord, const math::d3 &domain_size) {
    auto container = Registry::instance->moleculeContainer();

    const math::ul3 cell_dims = container->getCells().dims();
    math::d3 is_bound {0, 0, 0};
    if (cell_coord.x() == 1) is_bound.x() = 1;
    if (cell_coord.x() == cell_dims.x() - 2) is_bound.x() = -1;
    if (cell_coord.y() == 1) is_bound.y() = 1;
    if (cell_coord.y() == cell_dims.y() - 2) is_bound.y() = -1;
    if (cell_coord.z() == 1) is_bound.z() = 1;
    if (cell_coord.z() == cell_dims.z() - 2) is_bound.z() = -1;

    uint64_t idx = 0;
    for (int shiftZ = 0; shiftZ <= ((is_bound.z() != 0) ? 1 : 0); shiftZ++) {
        for (int shiftY = 0; shiftY <= ((is_bound.y() != 0) ? 1 : 0); shiftY++) {
            for (int shiftX = 0; shiftX <= ((is_bound.x() != 0) ? 1 : 0); shiftX++) {
                if (shiftX == 0 && shiftY == 0 && shiftZ == 0) continue;
                const math::i3 shiftVec {shiftX, shiftY, shiftZ};
                auto& [halo_cell_ref, halo_shift] = molecule.getCopies()[idx++];
                if (halo_shift != shiftVec) throw std::runtime_error("Something went wrong!");

                Cell& halo_cell = halo_cell_ref.get();
                Molecule& halo_molecule = *halo_cell.findBy(molecule.ID());
                const math::d3 halo_offset = domain_size * is_bound * shiftVec;
                halo_molecule.copyParentLocation(halo_offset);
            }
        }
    }
}

void Boundary::loopOverBoundary(std::function<void(Cell&, const math::ul3&)> fun) {
    auto container = Registry::instance->moleculeContainer();
    auto& cells = container->getCells();
    const math::ul3 cell_dims = container->getCells().dims();

    for (uint64_t z = 1; z < cell_dims.z()-1; z++) {
        for (uint64_t y = 1; y < cell_dims.y()-1; y++) {
            for (uint64_t x = 1; x < cell_dims.x()-1; x++) {
                if ((x != 1 && x != cell_dims.x()-2) &&
                    (y != 1 && y != cell_dims.y()-2) &&
                    (z != 1 && z != cell_dims.z()-2)) continue;

                Cell& cell = cells[x, y, z];
                fun(cell, {x, y, z});
            }
        }
    }
}
