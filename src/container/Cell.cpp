//
// Created by alex on 8/1/24.
//

#include "Cell.h"
#include "molecule/Molecule.h"
#include "SOA.h"
#include "math/Geometry.h"
#include "Registry.h"

Cell Cell::INVALID = Cell({1, 1, 1}, {0, 0, 0});

Cell::Cell(const math::d3 &low, const math::d3 &high) : m_low(low), m_high(high) { }

void Cell::setBounds(const math::d3 &low, const math::d3 &high) {
    m_low = low;
    m_high = high;
}

void Cell::addMolecule(const Molecule &molecule) {
    m_data.push_back(molecule);
    m_data.back().setCell(*this);
}

void Cell::constructSOA() {
    if (!Registry::instance->configuration()->enableSOA) return;
    // Allocate memory
    uint64_t num_sites = 0;
    for (Molecule& molecule : m_data) num_sites += molecule.getSites().size();
    m_soa.resize(num_sites);

    // write to soa
    uint64_t idx = 0;
    for (Molecule& molecule : m_data) {
        for (Site& site : molecule.getSites()) {
            m_soa.epsilon()[idx] = site.getEpsilon();
            m_soa.sigma()[idx] = site.getSigma();
            m_soa.mass()[idx] = site.getMass();
            m_soa.id()[idx] = molecule.ID();
            m_soa.r()[idx] = site.r_arr();
            m_soa.f()[idx] = site.f_arr();
            m_soa.fold()[idx] = site.fold_arr();
            m_soa.v()[idx] = site.v_arr();
            idx++;
        }
    }
    m_valid_soa = true;
}

void Cell::writeSOA2AOS() {
    if (!Registry::instance->configuration()->enableSOA) return;

    // the order should not have changed from writing to SOA
    uint64_t idx = 0;
    for (Molecule& molecule : m_data) {
        for (Site& site : molecule.getSites()) {
            site.r_arr() = m_soa.r()[idx];
            site.f_arr() = m_soa.f()[idx];
            site.fold_arr() = m_soa.fold()[idx];
            site.v_arr() = m_soa.v()[idx];
            idx++;
        }
    }
}

bool Cell::insideBounds(const math::d3 &point) {
    return math::pointInBox(point, m_low, m_high);
}

void Cell::invalidateSOA() {
    m_valid_soa = false;
}

std::vector<Molecule>::iterator Cell::findBy(uint64_t id) {
    for (auto it = m_data.begin(); it != m_data.end(); ++it) {
        if (it->ID() == id) return m_data.begin() + std::distance(m_data.begin(), it);
    }
    return m_data.end();
}

std::vector<Molecule>::iterator Cell::removeMolecule(uint64_t id) {
    auto it = findBy(id);
    if (it == m_data.end()) throw std::runtime_error("Molecule was not in this cell!");

    it->setCell(Cell::INVALID);
    return m_data.erase(it);
}
