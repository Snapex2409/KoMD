//
// Created by alex on 7/31/24.
//

#include "Molecule.h"
#include "container/Cell.h"

Molecule Molecule::INVALID = Molecule();

Molecule::id_t Molecule::NEXT_ID = 0;

Molecule::Molecule() : m_sites(), m_id(NEXT_ID++), m_cid() { }

Molecule::Molecule(const Molecule &copy) {
    m_sites = copy.m_sites;
    m_id = copy.m_id;
    m_cid = copy.m_cid;
}

Molecule::Molecule(Molecule &&move) noexcept {
    m_sites = std::move(move.m_sites);
    m_id = move.m_id;
    m_cid = move.m_cid;
}

Molecule& Molecule::operator=(const Molecule &copy) {
    m_sites = copy.m_sites;
    m_id = copy.m_id;
    m_cid = copy.m_cid;
    return *this;
}

Molecule& Molecule::operator=(Molecule &&move) noexcept {
    m_sites = std::move(move.m_sites);
    m_id = move.m_id;
    m_cid = move.m_cid;
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

void Molecule::moveBy(const math::d3 &offset) {
    for (Site& site : m_sites) {
        site.r_arr() += offset;
    }
}

void Molecule::moveCoMTo(const math::d3 &position) {
    const math::d3 delta = position - getCenterOfMass();
    moveBy(delta);
}
