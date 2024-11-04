//
// Created by alex on 11/3/24.
//

#include "TripleList.h"
#include "Registry.h"

TripleList::TripleList() {
    auto config = Registry::instance->configuration();
    m_update_freq = config->triple_steps;
    m_update_counter = 0;

    m_num_triplets = 0;
}

void TripleList::init(int num_sites) {
    m_num_triplets = 0;
    m_triplets = KW::nvec_t<int, 3>("TripleList - Triplets", 1);
    m_triplets_offsets = KW::nvec_t<math::d3, 3>("TripleList - Offsets", 1);
}

void TripleList::resize(uint64_t num_target_triplets) {
    Kokkos::resize(m_triplets, num_target_triplets);
    Kokkos::resize(m_triplets_offsets, num_target_triplets);
    m_num_triplets = num_target_triplets;
}

void TripleList::addTriplet(int idx0, int idx1, int idx2,
                            const math::d3 &offset0,
                            const math::d3 &offset1,
                            const math::d3 &offset2) {
    uint64_t updated_pos = Kokkos::atomic_add_fetch(&m_num_triplets, 1);
    const uint64_t local_pos = updated_pos - 1;

    m_triplets(local_pos, 0) = idx0;
    m_triplets(local_pos, 1) = idx1;
    m_triplets(local_pos, 2) = idx2;
    m_triplets_offsets(local_pos, 0) = offset0;
    m_triplets_offsets(local_pos, 1) = offset1;
    m_triplets_offsets(local_pos, 2) = offset2;
}
