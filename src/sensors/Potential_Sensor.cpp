//
// Created by alex on 8/8/24.
//

#include "Potential_Sensor.h"
#include "Registry.h"
#include "IO/Logging.h"

#include <utility>
#include <sstream>
#include <fstream>

Potential_Sensor::Potential_Sensor(std::string name, uint64_t bins) :
    Sensor(std::move(name)), p_bins(bins), p_data_u(bins, 0.0), p_data_f(bins, 0.0),
    p_count_u(bins, 0.0), p_count_f(bins, 0.0),
    m_sigma(Registry::instance->configuration()->sigma) { }

void Potential_Sensor::measure() {
    ForceFunctor::operator()();
}

void Potential_Sensor::write(uint64_t simstep) {
    std::stringstream file_name;
    file_name << "pot_" << name() << "_" << simstep << ".txt";
    std::ofstream file(file_name.str());
    if (!file.is_open()) {
        Log::io->error() << "Could not write file: " << file_name.str() << std::endl;
        return;
    }

    double r = 0;
    const double max_size = 3.0 * m_sigma;
    const double bin_width = max_size / static_cast<double>(p_bins);
    for (uint64_t idx = 0; idx < p_bins; idx++) {
        file << r << " " << p_data_f[idx] / (double) p_count_f[idx] << " " << p_data_u[idx] / (double) p_count_u[idx] << "\n";
        r += bin_width;
    }

    file.close();
}

uint64_t Potential_Sensor::get_bin(double r) {
    if (r < 0) throw std::runtime_error("negative r not allowed");

    const double max_size = 3.0 * m_sigma;
    const double bin_width = max_size / static_cast<double>(p_bins);
    auto bin = static_cast<uint64_t>(r / bin_width);
    return std::clamp(bin, 0UL, p_bins-1);
}
