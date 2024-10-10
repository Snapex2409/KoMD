//
// Created by alex on 10/10/24.
//

#ifndef KOMD_PAIRLIST_H
#define KOMD_PAIRLIST_H

#include "util/Kokkos_Wrapper.h"

class PairList {
public:
    /// allocates memory for all potential interactions, s.t. no reallocations will be needed
    void init(int num_sites);

    /// clears buffers by setting size counter to 0
    void reset() { m_num_pairs = 0; }

    /// insert one pair into the buffer
    void addPair(int idx0, int idx1, const math::d3& offset0 = {0, 0, 0}, const math::d3& offset1 = {0, 0, 0});

    KW::nvec_t<int, 2>& getPairs() { return m_pairs; }
    KW::nvec_t<math::d3, 2>& getOffsets() { return m_pair_offsets; }
    [[nodiscard]] uint64_t size() const { return m_num_pairs; }
private:
    /// index pairs
    KW::nvec_t<int, 2> m_pairs;
    /// offsets for pairs
    KW::nvec_t<math::d3, 2> m_pair_offsets;
    /// total number of active pairs in buffer
    uint64_t m_num_pairs;
};


#endif //KOMD_PAIRLIST_H
