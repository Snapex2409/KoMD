//
// Created by alex on 8/8/24.
//

#include "CheckpointIO.h"

#include "Registry.h"
#include "Logging.h"

#include <sstream>
#include <fstream>

void CheckpointIO::writeCheckpoint(uint64_t simstep) {
    auto container = Registry::instance->moleculeContainer();
    container->writeSOA2AOS();

    std::stringstream file_name;
    file_name << "ckpt_" << simstep << ".ps";
    std::ofstream file(file_name.str());
    if (!file.is_open()) {
        Log::io->error() << "Could not write checkpoint file: " << file_name.str() << std::endl;
        return;
    }

    auto& cells = container->getCells();
    const math::ul3 cell_dims = cells.dims();
    for (uint64_t z = 1; z < cell_dims.z()-1; z++) {
        for (uint64_t y = 1; y < cell_dims.y()-1; y++) {
            for (uint64_t x = 1; x < cell_dims.x()-1; x++) {
                Cell& cell = cells[x, y, z];
                for (Molecule& molecule : cell.molecules()) {
                    for (Site& site : molecule.getSites()) {
                        file << molecule.ID() << "\t";
                        file << site.r_arr().x() << " " << site.r_arr().y() << " " << site.r_arr().z() << "\t";
                        file << site.v_arr().x() << " " << site.v_arr().y() << " " << site.v_arr().z() << "\t";
                        file << site.getMass()   << " " << site.getEpsilon()<< " " << site.getSigma()  << "\n";
                    }
                }
            }
        }
    }

    file.close();
}

void CheckpointIO::loadCheckpoint() {
    auto config = Registry::instance->configuration();
    auto container = Registry::instance->moleculeContainer();
    const std::string& file_name = config->checkpoint_file;

    std::ifstream file(file_name);
    if (!file.is_open()) {
        Log::io->error() << "Could not open checkpoint file: " << file_name << std::endl;
        return;
    }

    uint64_t id;
    math::d3 r, v;
    double mass, eps, sig;
    bool read_line = false;
    Molecule molecule;
    while (!file.eof() && file.good()) {
        file >> id >> r.x() >> r.y() >> r.z() >> v.x() >> v.y() >> v.z() >> mass >> eps >> sig;

        if (!read_line) {
            molecule.setID(id);
            read_line = true;
        }

        if (id != molecule.ID()) {
            container->addMolecule(molecule);
            molecule = Molecule();
            molecule.setID(id);
        }
        molecule.addSite(eps, sig, mass, r, v);
    }

    // add last molecule
    if (read_line) container->addMolecule(molecule);

    file.close();
}
