//
// Created by alex on 8/22/24.
//

#include "Component.h"

Component::Component() : m_component_id(), m_sites() { }

void Component::addSite(double epsilon, double sigma, double mass, const math::d3 &r) {
    m_sites.emplace_back(epsilon, sigma, mass, r);
}

Molecule::sites_t &Component::getSites() {
    return m_sites;
}

