//
// Created by alex on 7/31/24.
//

#ifndef KOMD_MOLECULE_H
#define KOMD_MOLECULE_H

#include <vector>
#include <functional>

#include "math/Array.h"
#include "util/defaults.h"
#include "Site.h"

class Cell;

/**
 * Container for individual particles.
 * */
class Molecule {
public:
    using sites_t = std::vector<Site>;
    using id_t = uint32_t;

    Molecule();
    Molecule(const Molecule& copy);
    Molecule(Molecule&& move) noexcept;
    Molecule& operator=(const Molecule& copy);
    Molecule& operator=(Molecule&& move) noexcept;

    /**
     * Adds a new site at position r.
     * */
    void addSite(double epsilon = Defaults::epsilon, double sigma = Defaults::sigma, double mass = Defaults::mass, const math::d3& r = Defaults::r, const math::d3& v = Defaults::v);

    /**
     * Gets molecule id
     * */
    [[nodiscard]] id_t ID() const { return m_id; }

    /**
     * Sets the molecule id
     * */
    void setID(id_t id) { m_id = id; }

    /**
     * Gets molecule cid
     * */
    [[nodiscard]] id_t CID() const { return m_cid; }

    /**
     * Sets the molecule cid
     * */
    void setCID(id_t cid) { m_cid = cid; }

    /**
     * Gets reference to this indices sites
     * */
    sites_t& getSites();

    /**
     * Returns the center of mass position, following:\n
     * R = (sum m_i * r_i) / (sum m_i)
     * */
    [[nodiscard]] math::d3 getCenterOfMass() const;

    /**
     * Moves all sites of this molecule by offset\n
     * !! Not to offset !!
     * */
    void moveBy(const math::d3& offset);

    /**
     * Moves Center of Mass to the provided position\n
     * */
    void moveCoMTo(const math::d3& position);

    static Molecule INVALID;
private:
    /// Collection of all sites
    sites_t m_sites;
    /// Globally unique id
    id_t m_id;
    /// Component id
    id_t m_cid;
    /// Used for id assignment
    static id_t NEXT_ID;
};


#endif //KOMD_MOLECULE_H
