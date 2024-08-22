//
// Created by alex on 8/22/24.
//

#ifndef KOMD_COMPONENT_H
#define KOMD_COMPONENT_H

#include "Site.h"
#include "Molecule.h"
#include "util/defaults.h"
/**
 * Acts as a "description" as how to construct molecules of this kind
 * */
class Component {
public:
    Component();

    /**
     * Sets the id
     * */
    void setID(id_t id) { m_component_id = id; }

    /**
     * Adds a new site at position r.
     * */
    void addSite(double epsilon = Defaults::epsilon, double sigma = Defaults::sigma, double mass = Defaults::mass, const math::d3& r = Defaults::r);

    /**
     * Gets component id
     * */
    [[nodiscard]] id_t ID() const { return m_component_id; }

    /**
     * Gets reference to this components sites
     * */
    Molecule::sites_t& getSites();
private:
    /// container of all sites
    Molecule::sites_t m_sites;
    /// unique identifier of this component type
    Molecule::id_t m_component_id;
};


#endif //KOMD_COMPONENT_H
