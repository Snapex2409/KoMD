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
        Sensor(std::move(name)), p_bins(bins), p_data_u("POT_u", bins), p_data_f("POT_u", bins),
        p_count_u("POT_u", bins), p_count_f("POT_u", bins),
        p_max_sigma(Registry::instance->configuration()->max_sigma) {
    p_data_u_scatter = SOA::vec_scatter_t<double>(p_data_u);
    p_data_f_scatter = SOA::vec_scatter_t<double>(p_data_f);
    p_count_u_scatter = SOA::vec_scatter_t<uint64_t>(p_count_u);
    p_count_f_scatter = SOA::vec_scatter_t<uint64_t>(p_count_f);
    p_run_contribution = false;
}

void Potential_Sensor::measure() {
    ForceFunctor::operator()();
    Kokkos::Experimental::contribute(p_data_f, p_data_f_scatter);
    Kokkos::Experimental::contribute(p_data_u, p_data_u_scatter);
    Kokkos::Experimental::contribute(p_count_f, p_count_f_scatter);
    Kokkos::Experimental::contribute(p_count_u, p_count_u_scatter);
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
    const double max_size = 3.0 * p_max_sigma;
    const double bin_width = max_size / static_cast<double>(p_bins);
    for (uint64_t idx = 0; idx < p_bins; idx++) {
        file << r << " " << p_data_f[idx] / (double) p_count_f[idx] << " " << p_data_u[idx] / (double) p_count_u[idx] << "\n";
        r += bin_width;
    }

    file.close();
}
