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
m_bins(), m_samples(0) {
    const uint64_t num_bins = m_max_r / m_delta_r;
    m_bins = KW::vec_t<double>("RDF", num_bins);
    m_com = KW::vec_t<math::d3>("RDF CoM", Registry::instance->moleculeContainer()->getNumMolecules());
}

void RDFSensor::measure() {
    auto config = Registry::instance->configuration();
    auto container = Registry::instance->moleculeContainer();
    const auto count = container->getNumMolecules();
    container->getCenterOfMassPositions(m_com);
    const math::d3 dom_size = config->domainHigh - config->domainLow;

    Kokkos::parallel_for("RDF", Kokkos::MDRangePolicy({0, 0}, {count, count}), RDF_Kernel(m_com, m_bins, m_max_r, m_delta_r, static_cast<uint64_t>(m_max_r / m_delta_r), dom_size));
    Kokkos::fence("RDF - fence");

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
    const double n = Registry::instance->moleculeContainer()->getNumMolecules();
    const double n_total = n * (n-1);
    double r = m_delta_r / 2.0;
    const double dr_half = r;
    for (uint64_t idx = 0; idx < m_bins.size(); idx++) {
        const double n_bin = m_bins[idx] / m_samples;
        const double v_bin = 4.0 / 3.0 * M_PI * (std::pow(r+dr_half, 3) - std::pow(r-dr_half, 3));
        const double g_r = (n_bin * v_total) / (v_bin * n_total);
        file << r << "\t" << g_r << "\t" << n_bin << "\t" << v_bin << "\t" << n_total << "\t" << v_total << "\t" << (n_bin * v_total) << "\t" << (v_bin * n_total) << "\n";
        r += m_delta_r;
    }
    file.close();
}

void RDFSensor::RDF_Kernel::operator()(int idx_0, int idx_1) const {
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
