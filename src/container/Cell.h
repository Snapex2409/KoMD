//
// Created by alex on 8/1/24.
//

#ifndef KOMD_CELL_H
#define KOMD_CELL_H

#include "math/Array.h"
#include "util/Kokkos_Wrapper.h"

class Cell {
public:
    Cell() = default;
    Cell(const math::d3& low, const math::d3& high, const math::ul3& coord);
    void setBounds(const math::d3& low, const math::d3& high);
    void setCoords(const math::ul3& coord);

    /**
     * Checks if the point is within the bounds of this cell
     * */
    bool insideBounds(const math::d3& point);

    [[nodiscard]] const math::d3& low() const { return m_low; }
    [[nodiscard]] const math::d3& high() const { return m_high; }
    [[nodiscard]] const math::ul3& coord() const { return m_coord; }
    [[nodiscard]] KW::vec_t<uint64_t>& indices() { return m_data; }
    void addIndex(uint64_t idx);
    [[nodiscard]] uint64_t getNumIndices() { return m_num_indices; }
    void resetIndices();

    /**
     * Creates Kokkos buffers that stores buffers into global SOA
     * */
    void createIndexBuffers(uint64_t size);

    static Cell INVALID;
private:
    /// lower corner of cell
    math::d3 m_low;
    /// upper corner of cell
    math::d3 m_high;
    /// coordinate of cell
    math::ul3 m_coord;
    /// buffer of indices of indices inside active molecule container
    KW::vec_t<uint64_t> m_data;
    /// number of inserted indices
    uint64_t m_num_indices;
};


#endif //KOMD_CELL_H
