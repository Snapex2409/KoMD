//
// Created by alex on 8/1/24.
//

#include "Cell.h"
#include "molecule/Molecule.h"
#include "SOA.h"
#include "math/Geometry.h"
#include "Registry.h"

Cell Cell::INVALID = Cell({1, 1, 1}, {0, 0, 0}, {0, 0, 0});

Cell::Cell(const math::d3 &low, const math::d3 &high, const math::ul3& coord) : m_low(low), m_high(high), m_coord(coord) { }

void Cell::setBounds(const math::d3 &low, const math::d3 &high) {
    m_low = low;
    m_high = high;
}

void Cell::setCoords(const math::ul3 &coord) {
    m_coord = coord;
}

bool Cell::insideBounds(const math::d3 &point) {
    return math::pointInBox(point, m_low, m_high);
}

void Cell::createIndexBuffers(const uint64_t size) {
    m_data = KW::vec_t<uint64_t>("Cell Index Buffer", size);
}

void Cell::addIndex(uint64_t idx) {
    m_data[m_num_indices++] = idx;
}

void Cell::resetIndices() {
    m_num_indices = 0;
}
