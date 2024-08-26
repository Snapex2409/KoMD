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

    bool begin = true;
    for (auto it = container->iterator(MoleculeContainer::SITE, MoleculeContainer::DOMAIN); it->isValid(); ++(*it)) {
        if (begin) begin = false;
        else file << "\n";

        file << it->ID() << "\t";
        file << it->r().x() << " " << it->r().y()  << " " << it->r().z()  << "\t";
        file << it->v().x() << " " << it->v().y()  << " " << it->v().z()  << "\t";
        file << it->mass()  << " " << it->epsilon()<< " " << it->sigma();
    }

    file.close();
}

void CheckpointIO::loadCheckpoint() {
    auto config = Registry::instance->configuration();
    auto container = Registry::instance->moleculeContainer();

    auto& checkpoint_files = config->checkpoint_files;
    for (const auto& [path, offset] : checkpoint_files) {
        std::ifstream file(path);
        if (!file.is_open()) {
            Log::io->error() << "Could not open checkpoint file: " << path << std::endl;
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
            molecule.addSite(eps, sig, mass, r + offset, v);
        }

        // add last molecule
        if (read_line) container->addMolecule(molecule);

        file.close();
    }
}
