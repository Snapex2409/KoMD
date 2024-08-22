//
// Created by alex on 8/11/24.
//

#include "RDFSensor.h"
#include "Registry.h"
#include "IO/Logging.h"

#include <utility>
#include <sstream>
#include <fstream>

RDFSensor::RDFSensor() : Sensor("RDF"),
m_max_r(Registry::instance->configuration()->sensor_rdf_max),
m_delta_r(Registry::instance->configuration()->sensor_rdf_dr),
m_rho_0(Registry::instance->configuration()->density),
m_bins(), m_bins_scatter(), m_samples(0) {
    const uint64_t num_bins = m_max_r / m_delta_r;
    m_bins = SOA::vec_t<double>("RDF", num_bins);
    m_bins_scatter = SOA::vec_scatter_t<double>(m_bins);
}

void RDFSensor::measure() {
    auto container = Registry::instance->moleculeContainer();

    const auto count = container->getNumMolecules();
    SOA::vec_t<math::d3> positions("RDF CoM", count);
    container->getCenterOfMassPositions(positions);

    Kokkos::parallel_for("RDF", Kokkos::MDRangePolicy({0, 0}, {count, count}), RDF_Kernel(positions, m_bins_scatter, m_max_r, m_delta_r, static_cast<uint64_t>(m_max_r / m_delta_r)));
    Kokkos::fence("RDF - fence");
    Kokkos::Experimental::contribute(m_bins, m_bins_scatter);

    m_samples+=1;
}

void RDFSensor::write(uint64_t simstep) {
    std::stringstream file_name;
    file_name << "rdf_" << simstep << ".txt";
    std::ofstream file(file_name.str());
    if (!file.is_open()) {
        Log::io->error() << "Could not write file: " << file_name.str() << std::endl;
        return;
    }

    auto config = Registry::instance->configuration();
    const double v_total = (config->domainHigh - config->domainLow).product();
    const double n = m_rho_0 * v_total;
    const double n_total = 0.5 * n * (n-1);
    double r = m_delta_r / 2.0;
    for (uint64_t idx = 0; idx < m_bins.size(); idx++) {
        const double n_bin = m_bins[idx] / m_samples;
        const double v_bin = 4.0 * M_PI * m_delta_r * r * r;
        const double g_r = (n_bin * v_total) / (v_bin * n_total);
        file << r << " " << g_r << "\n";
        r += m_delta_r;
    }

    file.close();
}

void RDFSensor::RDF_Kernel::operator()(int idx_0, int idx_1) const {
    auto bin_access = bins_scatter.access();
    const math::d3 r0 = positions[idx_0];
    const math::d3 r1 = positions[idx_1];
    const double r = (r0 - r1).L2();
    if (r > max_r) return;

    const uint64_t bin = getBin(r, delta_r, bins);
    bin_access(bin) += 1;
}
