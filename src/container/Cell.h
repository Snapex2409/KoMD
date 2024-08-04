//
// Created by alex on 8/1/24.
//

#ifndef KOMD_CELL_H
#define KOMD_CELL_H

#include "math/Array.h"
#include "SOA.h"
#include <vector>
#include <memory>

class Molecule;

class Cell {
public:
    Cell() = default;
    Cell(const math::d3& low, const math::d3& high);
    void setBounds(const math::d3& low, const math::d3& high);
    void addMolecule(const Molecule& molecule);
    std::vector<Molecule>::iterator removeMolecule(uint64_t id);
    void constructSOA();
    void writeSOA2AOS();
    void invalidateSOA();

    /**
     * Checks if the point is within the bounds of this cell
     * */
    bool insideBounds(const math::d3& point);

    /**
     * Searches for the location of a molecule defined by its id.\n
     * Location expressed in terms of iterator of the local molecule vector.\n
     * If molecule was not found, then end iterator is returned.
     * */
    std::vector<Molecule>::iterator findBy(uint64_t id);

    [[nodiscard]] const math::d3& low() const { return m_low; }
    [[nodiscard]] const math::d3& high() const { return m_high; }
    [[nodiscard]] SOA& soa() { return m_soa; }
    [[nodiscard]] std::vector<Molecule>& molecules() { return m_data; }
    [[nodiscard]] bool validSOA() { return m_valid_soa; }

    static Cell INVALID;
private:

    math::d3 m_low;
    math::d3 m_high;
    std::vector<Molecule> m_data;
    SOA m_soa;
    bool m_valid_soa = false;
};


#endif //KOMD_CELL_H
