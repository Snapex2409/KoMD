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
    m_bins.resize(num_bins, 0.0);
}

void RDFSensor::measure() {
    auto container = Registry::instance->moleculeContainer();

    for (auto it0 = container->iterator(MoleculeContainer::MOLECULE, MoleculeContainer::DOMAIN); it0.isValid(); ++it0) {
        auto it1 = it0;
        ++it1;
        for (; it1.isValid(); ++it1) {
            const math::d3 r0 = it0.molecule().getCenterOfMass();
            const math::d3 r1 = it1.molecule().getCenterOfMass();
            const double r = (r0 - r1).L2();
            if (r > m_max_r) continue;

            const uint64_t bin = getBin(r);
            m_bins[bin] += 2;
        }
    }

    m_samples+=1;
}

void RDFSensor::write(uint64_t simstep) {
    const double fac = 4.0 * M_PI * m_delta_r * m_rho_0;

    std::stringstream file_name;
    file_name << "rdf_" << simstep << ".txt";
    std::ofstream file(file_name.str());
    if (!file.is_open()) {
        Log::io->error() << "Could not write file: " << file_name.str() << std::endl;
        return;
    }

    double r = m_delta_r / 2.0;
    for (uint64_t idx = 0; idx < m_bins.size(); idx++) {
        const double g_r = (m_bins[idx] / m_samples) / (r * r * fac);
        file << r << " " << g_r << "\n";
        r += m_delta_r;
    }

    file.close();
}

uint64_t RDFSensor::getBin(double r) const {
    if (r < 0) throw std::runtime_error("negative r not allowed");

    auto bin = static_cast<uint64_t>(r / m_delta_r);
    return std::clamp(bin, 0UL, m_bins.size()-1);
}
