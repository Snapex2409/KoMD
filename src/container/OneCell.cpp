//
// Created by alex on 8/20/24.
//

#include "OneCell.h"
#include "Registry.h"
#include "math/Geometry.h"

OneCell::OneCell() : MoleculeContainer(), m_data() {
    auto config = Registry::instance->configuration();
    m_data.setBounds(config->domainLow, config->domainHigh);
}

void OneCell::init() {
    uint64_t num_sites = 0;
    for (auto& molecule : p_molecules) num_sites += molecule.getSites().size();
    p_soa.createBuffers(num_sites);
    constructSOAs();
    m_data.createIndexBuffers(num_sites);
    for (uint64_t idx = 0; idx < num_sites; idx++) m_data.addIndex(idx);

    p_com = KW::vec_t<math::d3>("Center of Masses", p_soa.size());
}

void OneCell::updateContainer() {
    updateCOM();

    // apply periodic boundary kernel
    Kokkos::parallel_for("Periodic Bound", p_soa.size(), OneCell::Periodic_Kernel(p_com, p_soa.r(), m_data.low(), m_data.high(), m_data.high() - m_data.low()));
    // apply force reset kernel
    Kokkos::parallel_for("Force Reset", p_soa.size(), OneCell::FReset_Kernel(p_soa.f()));

    Kokkos::fence("OneCell - update");

    // update com vector
    updateCOM();
}

std::unique_ptr<MoleculeContainer::CellIterator> OneCell::iteratorCell() {
    return std::make_unique<OCCellIterator>(m_data);
}

std::unique_ptr<MoleculeContainer::CellPairIterator> OneCell::iteratorC08() {
    return std::make_unique<OCC08Iterator>(m_data);
}

void OneCell::FReset_Kernel::operator()(int idx) const {
    f[idx] = math::d3 {0, 0, 0};
}

void OneCell::Periodic_Kernel::operator()(int idx) const {
    const math::d3 pos = com[idx];

    math::d3 offset {0, 0, 0};
    for (int dim = 0; dim < 3; dim++) {
        if (pos[dim] < low[dim]) offset[dim] = domain_size[dim];
        if (pos[dim] > high[dim]) offset[dim] = -domain_size[dim];
    }

    r[idx] += offset;
}
