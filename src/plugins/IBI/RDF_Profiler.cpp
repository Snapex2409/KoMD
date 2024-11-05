//
// Created by alex on 11/4/24.
//

#include "RDF_Profiler.h"
#include "Registry.h"

RDF_Profiler::RDF_Profiler() : m_samples(0) {
    m_max_r = Registry::instance->configuration()->cutoff;
    const uint64_t num_bins = Registry::instance->configuration()->IBI_bins;
    m_delta_r = m_max_r / static_cast<double>(num_bins);
    m_bins = KW::vec_t<double>("IBI RDF", num_bins);
    m_com = KW::vec_t<math::d3>("IBI RDF CoM", Registry::instance->moleculeContainer()->getNumMolecules());
    m_nodes = KW::vec_t<double>("IBI nodes", num_bins);
    for (int i = 0; i < num_bins; ++i) m_nodes[i] = (i + 0.5) * m_delta_r;
}

void RDF_Profiler::measure() {
    auto config = Registry::instance->configuration();
    auto container = Registry::instance->moleculeContainer();
    const auto count = container->getNumMolecules();
    container->getCenterOfMassPositions(m_com);
    const math::d3 dom_size = config->domainHigh - config->domainLow;

    Kokkos::parallel_for("IBI RDF", Kokkos::MDRangePolicy({0, 0}, {count, count}),
        IBI_RDF_Kernel(m_com, m_bins, m_max_r, m_delta_r, m_bins.size(), dom_size));
    Kokkos::fence("IBI RDF - fence");

    m_samples+=1;
}

void RDF_Profiler::reset() {
    for (int idx = 0; idx < m_nodes.size(); ++idx) m_bins[idx] = 0;
    m_samples = 0;
}

void RDF_Profiler::getRDF(KW::vec_t<double> &buffer) const {
    // copy current samples into buffer
    for (int idx = 0; idx < m_nodes.size(); ++idx) buffer[idx] = m_bins[idx];

    // normalize
    auto config = Registry::instance->configuration();
    const double N = static_cast<double>(Registry::instance->moleculeContainer()->getNumMolecules());
    const double V = (config->domainHigh - config->domainLow).product();

    for (int bin = 0; bin < m_nodes.size(); ++bin) {
        const auto  rmin = bin * m_delta_r;
        const auto  rmax =(bin + 1) * m_delta_r;
        const auto  rmin3 = rmin * rmin * rmin;
        const auto  rmax3 = rmax * rmax * rmax;
        const auto  binvol = (4.0 / 3.0) * M_PI * (rmax3 - rmin3);
        const auto  den = N * (N - 1.0) * binvol / V;
        buffer[bin] /= m_samples;
        buffer[bin] /= den;
    }
}

uint64_t RDF_Profiler::getBufferSize() const {
    return m_nodes.size();
}

void RDF_Profiler::getRNodes(KW::vec_t<double> &buffer) const {
    for (int idx = 0; idx < m_nodes.size(); ++idx) buffer[idx] = m_nodes[idx];
}

void RDF_Profiler::IBI_RDF_Kernel::operator()(int idx_0, int idx_1) const {
    if (idx_0 == idx_1) return;
    const math::d3 r0 = positions[idx_0];
    const math::d3 r1 = positions[idx_1];
    const math::d3 delta = r0 - r1;
    const math::d3 dr = delta - domain_size * math::round(delta / domain_size);
    const double r = dr.L2();
    if (r > max_r) return;

    const uint64_t bin = getBin(r, delta_r, num_bins);
    Kokkos::atomic_add(&(bins[bin]), 1);
}
