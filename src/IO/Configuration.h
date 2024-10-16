//
// Created by alex on 7/31/24.
//

#ifndef KOMD_CONFIGURATION_H
#define KOMD_CONFIGURATION_H

#include "util/defaults.h"
#include "math/Array.h"
#include <string>
#include <vector>

struct Configuration {
    double temperature = Defaults::temperature;
    double delta_t = Defaults::delta_t;
    uint64_t timesteps = Defaults::timesteps;
    double cutoff = Defaults::cutoff;
    double cell_size = Defaults::cell_size;
    double density = Defaults::density;
    double max_epsilon = std::numeric_limits<double>::min();
    double max_sigma = std::numeric_limits<double>::min();
    double stiffness_factor = Defaults::stiffness_factor;
    double limit_factor = Defaults::limit_factor;

    int pair_steps = Defaults::pair_steps;

    math::d3 domainLow = Defaults::domainLow;
    math::d3 domainHigh = Defaults::domainHigh;
    std::vector<std::pair<std::string, math::d3>> checkpoint_files;
    uint64_t write_freq = Defaults::write_freq;
    /// begin, end, cid
    std::vector<std::tuple<math::d3, math::d3, uint32_t>> phasespace_gen_regions;

    bool storeCheckpoint = Defaults::storeCheckpoint;
    bool enable_sensor_lj = Defaults::enable_sensor_lj;
    bool enable_sensor_fene = Defaults::enable_sensor_fene;
    bool enable_sensor_rdf = Defaults::enable_sensor_rdf;
    bool enable_sensor_disp = Defaults::enable_sensor_disp;

    uint64_t sensor_lj_bins = Defaults::sensor_lj_bins;
    uint64_t sensor_fene_bins = Defaults::sensor_fene_bins;
    bool sensor_temp_opt = Defaults::sensor_temp_opt;
    double sensor_rdf_max = Defaults::sensor_rdf_max;
    double sensor_rdf_dr = Defaults::sensor_rdf_dr;

    bool ADR_enable = false;
    math::d3 ADR_low = Defaults::ADR_low;
    math::d3 ADR_high = Defaults::ADR_high;
    math::d3 ADR_h_dim = Defaults::ADR_h_dim;

    bool enable_one_cell = false;
};


#endif //KOMD_CONFIGURATION_H
