//
// Created by alex on 7/31/24.
//

#include "Molecule.h"
#include "container/Cell.h"

Molecule Molecule::INVALID = Molecule();

Molecule::id_t Molecule::NEXT_ID = 0;

Molecule::Molecule() : m_sites(), m_id(NEXT_ID++), m_parent(Cell::INVALID), m_links(), m_cell(Cell::INVALID) {}

Molecule::Molecule(const Molecule &copy) : m_parent(copy.m_parent), m_cell(copy.m_cell) {
    m_sites = copy.m_sites;
    m_id = copy.m_id;
    m_links = copy.m_links;
}

Molecule::Molecule(Molecule &&move) noexcept : m_parent(move.m_parent), m_cell(move.m_cell) {
    m_sites = std::move(move.m_sites);
    m_id = move.m_id;
    m_links = std::move(move.m_links);
}

Molecule& Molecule::operator=(const Molecule &copy) {
    m_sites = copy.m_sites;
    m_id = copy.m_id;
    m_parent = copy.m_parent;
    m_links = copy.m_links;
    m_cell = copy.m_cell;
    return *this;
}

Molecule& Molecule::operator=(Molecule &&move) noexcept {
    m_sites = std::move(move.m_sites);
    m_id = move.m_id;
    m_parent = move.m_parent;
    m_links = std::move(move.m_links);
    m_cell = move.m_cell;
    return *this;
}

void Molecule::addSite(double epsilon, double sigma, double mass, const math::d3 &r, const math::d3 &v) {
    m_sites.emplace_back(epsilon, sigma, mass, r, v);
}

Molecule::sites_t &Molecule::getSites() {
    return m_sites;
}

math::d3 Molecule::getCenterOfMass() const {
    math::d3 result {0, 0, 0};
    double total_mass = 0;
    for (const Site& site : m_sites) {
        result += site.r_arr() * site.getMass();
        total_mass += site.getMass();
    }
    return result / total_mass;
}

void Molecule::registerCopy(Cell& cell, const math::i3& shift) {
    m_links.emplace_back(cell, shift);
}

std::vector<std::pair<std::reference_wrapper<Cell>, math::i3>> &Molecule::getCopies() {
    return m_links;
}

void Molecule::setParent(Cell& cell) {
    m_parent = std::reference_wrapper<Cell>(cell);
}

Cell &Molecule::getParent() {
    return m_parent.get();
}

void Molecule::setCell(Cell &cell) {
    m_cell = std::reference_wrapper<Cell>(cell);
}

Cell &Molecule::getCell() {
    return m_cell.get();
}

void Molecule::moveBy(const math::d3 &offset) {
    for (Site& site : m_sites) {
        site.r_arr() += offset;
    }
}

void Molecule::moveCoMTo(const math::d3 &position) {
    const math::d3 delta = position - getCenterOfMass();
    moveBy(delta);
}

void Molecule::copyParentLocation(const math::d3 &offset) {
    for (uint64_t idx = 0; idx < m_sites.size(); idx++) {
        Site& site = m_sites[idx];
        Site& p_site = m_parent.get().findBy(m_id)->m_sites[idx];
        site.r_arr() = p_site.r_arr() + offset;
    }
}
