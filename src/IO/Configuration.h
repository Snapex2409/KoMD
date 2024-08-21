//
// Created by alex on 7/31/24.
//

#ifndef KOMD_CONFIGURATION_H
#define KOMD_CONFIGURATION_H

#include "util/defaults.h"
#include "math/Array.h"
#include <string>

struct Configuration {
    double temperature = Defaults::temperature;
    double delta_t = Defaults::delta_t;
    double cutoff = Defaults::cutoff;
    double density = Defaults::density;
    double epsilon = Defaults::epsilon;
    double sigma = Defaults::sigma;
    double stiffness_factor = Defaults::stiffness_factor;
    double limit_factor = Defaults::limit_factor;
    math::d3 domainLow = Defaults::domainLow;
    math::d3 domainHigh = Defaults::domainHigh;
    bool loadCheckpoint = Defaults::loadCheckpoint;
    bool storeCheckpoint = Defaults::storeCheckpoint;
    bool enable_sensor_lj = Defaults::enable_sensor_lj;
    bool enable_sensor_fene = Defaults::enable_sensor_fene;
    bool enable_sensor_rdf = Defaults::enable_sensor_rdf;
    uint64_t timesteps = Defaults::timesteps;
    uint64_t write_freq = Defaults::write_freq;
    std::string checkpoint_file = Defaults::checkpoint_file;

    uint64_t sensor_lj_bins = Defaults::sensor_lj_bins;
    uint64_t sensor_fene_bins = Defaults::sensor_fene_bins;
    bool sensor_temp_opt = Defaults::sensor_temp_opt;
    double sensor_rdf_max = Defaults::sensor_rdf_max;
    double sensor_rdf_dr = Defaults::sensor_rdf_dr;

    bool enable_one_cell = false;
};


#endif //KOMD_CONFIGURATION_H
