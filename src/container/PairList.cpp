//
// Created by alex on 10/10/24.
//

#include "PairList.h"
#include "Registry.h"

PairList::PairList() {
    auto config = Registry::instance->configuration();
    m_update_freq = config->pair_steps;
    m_update_counter = 0;

    m_num_pairs = 0;
}

void PairList::init(int num_sites) {
    m_num_pairs = 0;
    m_pairs = KW::nvec_t<int, 2>("PairList - Pairs", 1);
    m_pair_offsets = KW::nvec_t<math::d3, 2>("PairList - Offsets", 1);
}

void PairList::resize(uint64_t num_target_pairs) {
    Kokkos::resize(m_pairs, num_target_pairs);
    Kokkos::resize(m_pair_offsets, num_target_pairs);
    m_num_pairs = num_target_pairs;
}

void PairList::addPair(int idx0, int idx1, const math::d3 &offset0, const math::d3 &offset1) {
    uint64_t updated_pos = Kokkos::atomic_add_fetch(&m_num_pairs, 1);
    const uint64_t local_pos = updated_pos - 1;

    m_pairs(local_pos, 0) = idx0;
    m_pairs(local_pos, 1) = idx1;
    m_pair_offsets(local_pos, 0) = offset0;
    m_pair_offsets(local_pos, 1) = offset1;
}
