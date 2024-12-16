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
    double energy_3b = Defaults::energy_3b;

    int pair_steps = Defaults::pair_steps;
    int triple_steps = Defaults::triple_steps;

    math::d3 domainLow = Defaults::domainLow;
    math::d3 domainHigh = Defaults::domainHigh;
    std::vector<std::pair<std::string, math::d3>> checkpoint_files;
    uint64_t write_freq = Defaults::write_freq;
    /// begin, end, cid
    std::vector<std::tuple<math::d3, math::d3, uint32_t>> phasespace_gen_regions;

    bool storeCheckpoint = Defaults::storeCheckpoint;
    bool enable_3b_approx = Defaults::enable_3b_approx;
    bool enable_3b = Defaults::enable_3b;
    bool enable_3b_direct = Defaults::enable_3b_direct;
    bool enable_sensor_lj = Defaults::enable_sensor_lj;
    bool enable_sensor_fene = Defaults::enable_sensor_fene;
    bool enable_sensor_rdf = Defaults::enable_sensor_rdf;
    bool enable_sensor_disp = Defaults::enable_sensor_disp;
    bool enable_sensor_pres = Defaults::enable_sensor_pres;
    bool enable_sensor_visc = Defaults::enable_sensor_visc;

    uint64_t sensor_lj_bins = Defaults::sensor_lj_bins;
    uint64_t sensor_fene_bins = Defaults::sensor_fene_bins;
    bool sensor_temp_opt = Defaults::sensor_temp_opt;
    double sensor_rdf_max = Defaults::sensor_rdf_max;
    double sensor_rdf_dr = Defaults::sensor_rdf_dr;
    int sensor_visc_window = Defaults::sensor_visc_window;
    bool sensor_visc_no_orig = Defaults::sensor_visc_no_orig;

    bool ADR_enable = false;
    math::d3 ADR_low = Defaults::ADR_low;
    math::d3 ADR_high = Defaults::ADR_high;
    math::d3 ADR_h_dim = Defaults::ADR_h_dim;

    bool IBI_enable = false;
    uint64_t IBI_bins = Defaults::IBI_bins;
    double IBI_alpha = Defaults::IBI_alpha;
    int IBI_steps_equil = Defaults::IBI_steps_equil;
    int IBI_steps_measure = Defaults::IBI_steps_measure;
    double IBI_conv_threshold = Defaults::IBI_conv_threshold;
    std::string IBI_conv_mode = Defaults::IBI_conv_mode;
    std::string IBI_conv_stop = Defaults::IBI_conv_stop;
    int IBI_conv_window = Defaults::IBI_conv_window;

    bool IBI_reload_enable = false;
    std::string IBI_reload_fpath = "";
    std::string IBI_reload_ppath = "";
    math::d3 IBI_exclusion_low = Defaults::IBI_exclusion_low;
    math::d3 IBI_exclusion_high = Defaults::IBI_exclusion_high;

    bool enable_one_cell = false;
};


#endif //KOMD_CONFIGURATION_H
