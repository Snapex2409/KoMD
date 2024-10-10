//
// Created by alex on 10/10/24.
//

#include "PairList.h"

void PairList::init(int num_sites) {
    m_num_pairs = 0;
    m_pairs = KW::nvec_t<int, 2>("PairList - Pairs", num_sites * num_sites * 8 * 8);
    m_pair_offsets = KW::nvec_t<math::d3, 2>("PairList - Offsets", num_sites * num_sites * 8 * 8);
}

void PairList::addPair(int idx0, int idx1, const math::d3 &offset0, const math::d3 &offset1) {
    m_pairs(m_num_pairs, 0) = idx0;
    m_pairs(m_num_pairs, 1) = idx1;
    m_pair_offsets(m_num_pairs, 0) = offset0;
    m_pair_offsets(m_num_pairs, 1) = offset1;
    m_num_pairs++;
}
