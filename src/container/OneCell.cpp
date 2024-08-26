//
// Created by alex on 8/20/24.
//

#include "OneCell.h"
#include "Registry.h"
#include "math/Geometry.h"

OneCell::OneCell() : MoleculeContainer(), m_data(), m_first_update(true) {
    auto config = Registry::instance->configuration();
    m_data.setBounds(config->domainLow, config->domainHigh);
}

void OneCell::addMolecule(const Molecule &molecule) {
    const math::d3 pos = molecule.getCenterOfMass();
    if (!math::pointInBox(pos, m_data.low(), m_data.high())) {
        throw std::runtime_error("Inserting Molecule at invalid position.");
    }

    m_data.molecules().push_back(molecule);
    p_molecule_count++;
}

void OneCell::updateContainer() {
    if (m_first_update) {
        m_data.soa().createBuffers();
        constructSOAs();
        m_first_update = false;
        m_com = KW::vec_t<math::d3>("Center of Masses", m_data.soa().size());
    }

    updateCOM();

    // apply periodic boundary kernel
    Kokkos::parallel_for("Periodic Bound", m_data.soa().size(), Periodic_Kernel(m_com, m_data.soa().r(), m_data.low(), m_data.high(), m_data.high() - m_data.low()));
    // apply force reset kernel
    Kokkos::parallel_for("Force Reset", m_data.soa().size(), FReset_Kernel(m_data.soa().f()));

    Kokkos::fence("OneCell - update");

    // update com vector
    updateCOM();
}

void OneCell::getCenterOfMassPositions(KW::vec_t<math::d3> &buffer) {
    uint64_t idx = 0;
    for (auto it = iterator(MOLECULE, DOMAIN); it->isValid(); ++(*it)) {
        buffer[idx] = it->molecule().getCenterOfMass();
        idx++;
    }
}

std::unique_ptr<MoleculeContainer::Iterator>
OneCell::iterator(const MoleculeContainer::IteratorType type, const MoleculeContainer::IteratorRegion region) {
    if (type == MOLECULE) writeSOA2AOS();
    if (region == DOMAIN_HALO) throw std::runtime_error("not supported");

    return std::make_unique<OCIterator>(type == MOLECULE, m_data);
}
std::unique_ptr<MoleculeContainer::CellIterator> OneCell::iteratorCell(const MoleculeContainer::IteratorRegion region) {
    return std::make_unique<OCCellIterator>(m_data);
}
std::unique_ptr<MoleculeContainer::CellPairIterator> OneCell::iteratorC08() {
    return std::make_unique<OCC08Iterator>(m_data);
}

void OneCell::writeSOA2AOS() { m_data.writeSOA2AOS(); }
void OneCell::constructSOAs() { m_data.constructSOA(); }
void OneCell::clearForces() { m_data.clearForces(); }
void OneCell::updateCOM() {
    uint64_t s_idx = 0;
    for (uint64_t idx = 0; idx < p_molecule_count; idx++) {
        const uint64_t num_sites = m_data.molecules()[idx].getSites().size();

        // compute com
        math::d3 com {0, 0, 0};
        double total_mass = 0;
        for (uint64_t s_counter = 0; s_counter < num_sites; s_counter++) {
            const double mass = m_data.soa().mass()[s_idx + s_counter];
            com += m_data.soa().r()[s_idx + s_counter] * mass;
            total_mass += mass;
        }
        com /= total_mass;

        // update com buffer
        for (uint64_t s_counter = 0; s_counter < num_sites; s_counter++) {
            m_com[s_idx++] = com;
        }
    }
}

OneCell::OCIterator::OCIterator(bool only_mol, Cell &cell) : m_only_molecule(only_mol), m_site_idx(0), m_molecule_idx(0), m_visited_sites_cell(0), m_cell(cell) {}

void OneCell::OCIterator::operator++() {
    auto& molecules = m_cell.molecules();
    // move to next site
    m_visited_sites_cell++;
    if (!m_only_molecule) {
        if (m_site_idx < molecules[m_molecule_idx].getSites().size() - 1) {
            m_site_idx++;
            return;
        }
    }

    // move to next molecule
    m_molecule_idx++;
    m_site_idx = 0;
}

bool OneCell::OCIterator::isValid() const { return m_molecule_idx < m_cell.molecules().size(); }

math::d3 &OneCell::OCIterator::f() { return m_cell.soa().f()[m_visited_sites_cell]; }
math::d3 &OneCell::OCIterator::r() { return m_cell.soa().r()[m_visited_sites_cell]; }
math::d3 &OneCell::OCIterator::v() { return m_cell.soa().v()[m_visited_sites_cell]; }
double OneCell::OCIterator::epsilon() { return m_cell.soa().epsilon()[m_visited_sites_cell]; }
double OneCell::OCIterator::sigma() { return m_cell.soa().sigma()[m_visited_sites_cell]; }
double OneCell::OCIterator::mass() { return m_cell.soa().mass()[m_visited_sites_cell]; }
uint64_t OneCell::OCIterator::ID() { return m_cell.soa().id()[m_visited_sites_cell]; }
Molecule &OneCell::OCIterator::molecule() { return m_cell.molecules()[m_molecule_idx]; }

void OneCell::Periodic_Kernel::operator()(int idx) const {
    const math::d3 pos = com[idx];

    math::d3 offset {0, 0, 0};
    for (int dim = 0; dim < 3; dim++) {
        if (pos[dim] < low[dim]) offset[dim] = domain_size[dim];
        if (pos[dim] > high[dim]) offset[dim] = -domain_size[dim];
    }

    r[idx] += offset;
}

void OneCell::FReset_Kernel::operator()(int idx) const {
    f[idx] = math::d3 {0, 0, 0};
}
