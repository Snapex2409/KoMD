//
// Created by alex on 7/31/24.
//

#include "MoleculeContainer.h"
#include "math/Geometry.h"
#include "Registry.h"

MoleculeContainer::MoleculeContainer() :  p_molecule_count(0), p_molecules(), p_soa(), p_pair_list(), p_triple_list() { }

void MoleculeContainer::addMolecule(const Molecule &molecule) {
    auto config = Registry::instance->configuration();
    const math::d3 pos = molecule.getCenterOfMass();
    if (!math::pointInBox(pos, config->domainLow, config->domainHigh)) {
        throw std::runtime_error("Inserting Molecule at invalid position.");
    }

    p_molecules.push_back(molecule);
    p_molecule_count++;
}

void MoleculeContainer::writeSOA2AOS() {
    uint64_t s_idx = 0;
    for (uint64_t m_idx = 0; m_idx < p_molecule_count; m_idx++) {
        Molecule& molecule = p_molecules[m_idx];
        for (Site& site : molecule.getSites()) {
            site.r_arr() = p_soa.r()[s_idx];
            site.f_arr() = p_soa.f()[s_idx];
            site.v_arr() = p_soa.v()[s_idx];
            s_idx++;
        }
    }
}

void MoleculeContainer::constructSOAs() {
    uint64_t s_idx = 0;
    for (uint64_t m_idx = 0; m_idx < p_molecule_count; m_idx++) {
        Molecule& molecule = p_molecules[m_idx];
        for (Site& site : molecule.getSites()) {
            p_soa.epsilon()[s_idx] = site.getEpsilon();
            p_soa.sigma()[s_idx] = site.getSigma();
            p_soa.mass()[s_idx] = site.getMass();
            p_soa.id()[s_idx] = molecule.ID();
            p_soa.idx()[s_idx] = m_idx;
            p_soa.r()[s_idx] = site.r_arr();
            p_soa.f()[s_idx] = site.f_arr();
            p_soa.v()[s_idx] = site.v_arr();
            s_idx++;
        }
    }
}

void MoleculeContainer::updateCOM() {
    uint64_t s_idx = 0;
    for (uint64_t idx = 0; idx < p_molecule_count; idx++) {
        const uint64_t num_sites = p_molecules[idx].getSites().size();

        // compute com
        math::d3 com {0, 0, 0};
        double total_mass = 0;
        for (uint64_t s_counter = 0; s_counter < num_sites; s_counter++) {
            const double mass = p_soa.mass()[s_idx + s_counter];
            com += p_soa.r()[s_idx + s_counter] * mass;
            total_mass += mass;
        }
        com /= total_mass;

        // update com buffer
        for (uint64_t s_counter = 0; s_counter < num_sites; s_counter++) {
            p_com[s_idx++] = com;
        }
    }
}

void MoleculeContainer::getCenterOfMassPositions(KW::vec_t<math::d3> &buffer) {
    updateCOM();

    uint64_t s_idx = 0;
    for (uint64_t m_idx = 0; m_idx < p_molecule_count; m_idx++) {
        Molecule& molecule = p_molecules[m_idx];
        buffer[m_idx] = p_com[s_idx];
        s_idx += molecule.getSites().size();
    }
}

//===========================================================================
// Iterator
//===========================================================================

MoleculeContainer::Iterator MoleculeContainer::iterator(MoleculeContainer::IteratorType type) {
    if (type == MOLECULE) writeSOA2AOS();
    return {type == MOLECULE, p_molecules, p_soa};
}

MoleculeContainer::Iterator::Iterator(bool only_mol, std::vector<Molecule>& molecules, SOA& soa) : m_only_molecule(only_mol), m_site_idx(0), m_molecule_idx(0), m_visited_sites(0), m_molecules(molecules), m_soa(soa) {}

void MoleculeContainer::Iterator::operator++() {
    // move to next site
    m_visited_sites++;
    if (!m_only_molecule) {
        if (m_site_idx < m_molecules[m_molecule_idx].getSites().size() - 1) {
            m_site_idx++;
            return;
        }
    }

    // move to next molecule
    m_molecule_idx++;
    m_site_idx = 0;
}

bool MoleculeContainer::Iterator::isValid() const { return m_molecule_idx < m_molecules.size(); }
math::d3 &MoleculeContainer::Iterator::f() { return m_soa.f()[m_visited_sites]; }
math::d3 &MoleculeContainer::Iterator::r() { return m_soa.r()[m_visited_sites]; }
math::d3 &MoleculeContainer::Iterator::v() { return m_soa.v()[m_visited_sites]; }
double MoleculeContainer::Iterator::epsilon() { return m_soa.epsilon()[m_visited_sites]; }
double MoleculeContainer::Iterator::sigma() { return m_soa.sigma()[m_visited_sites]; }
double MoleculeContainer::Iterator::mass() { return m_soa.mass()[m_visited_sites]; }
uint64_t MoleculeContainer::Iterator::ID() { return m_soa.id()[m_visited_sites]; }
Molecule &MoleculeContainer::Iterator::molecule() { return m_molecules[m_molecule_idx]; }
