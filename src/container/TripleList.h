//
// Created by alex on 11/3/24.
//

#ifndef KOMD_TRIPLELIST_H
#define KOMD_TRIPLELIST_H

#include "util/Kokkos_Wrapper.h"

class TripleList {
public:
    TripleList();
    /// allocates memory for all potential interactions, s.t. no reallocations will be needed
    void init(int num_sites);

    /// clears buffers by setting size counter to 0
    void reset() { m_num_triplets = 0; }

    /// reserves memory for num_target_triplets entries, all must be populated for correct execution
    void resize(uint64_t num_target_triplets);

    void step() { m_update_counter = (m_update_counter + 1) % m_update_freq; }

    bool requiresUpdate() const { return m_update_counter == 0; }

    /// insert one triplet into the buffer, method is MT-Safe
    void addTriplet(int idx0, int idx1, int idx2, const math::d3& offset0 = {0, 0, 0}, const math::d3& offset1 = {0, 0, 0}, const math::d3& offset2 = {0, 0, 0});

    KW::nvec_t<int, 3>& getTriplets() { return m_triplets; }
    KW::nvec_t<math::d3, 3>& getOffsets() { return m_triplets_offsets; }
    [[nodiscard]] uint64_t size() const { return m_num_triplets; }
private:
    /// index triplet
    KW::nvec_t<int, 3> m_triplets;
    /// offsets for triplets
    KW::nvec_t<math::d3, 3> m_triplets_offsets;
    /// total number of active triplets in buffer
    uint64_t m_num_triplets;
    /// update frequency
    int m_update_freq;
    /// update counter
    int m_update_counter;
};


#endif //KOMD_TRIPLELIST_H
