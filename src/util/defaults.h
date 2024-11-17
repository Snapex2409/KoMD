//
// Created by alex on 7/31/24.
//

#ifndef KOMD_DEFAULTS_H
#define KOMD_DEFAULTS_H

#include "math/Array.h"
#include <string>

namespace Defaults {
    static constexpr double sigma = 1.0;                // 1 Armstrong
    static constexpr double epsilon = 1.0;              // 1 Joule (from outside) (internally in u * A² / ps²)
    static constexpr double mass = 1.0;                 // 1 Dalton = 1u
    static constexpr double delta_t = 0.001;            // 1 ps, 0.001 = 1 fs
    static constexpr double stiffness_factor = 30;      // unit-less
    static constexpr double limit_factor = 20;          // unit-less
    static constexpr double energy_3b = 1;              // 1 Joule*nm**9 from outside (internally in A**9 * u * A² / ps²)
    static constexpr double density = 0.1;              // unit-less, number density
    static constexpr math::d3 r {0, 0, 0};
    static constexpr math::d3 v {0, 0, 0};
    static constexpr double temperature = 300;          // 300 Kelvin
    static constexpr math::d3 domainLow {0, 0, 0};      // in Armstrong
    static constexpr math::d3 domainHigh {10, 10, 10};  // in Armstrong
    static constexpr double cutoff = sigma * 2.5;       // in Armstrong
    static constexpr double cell_size = cutoff * 3.0 / 2.5;       // in Armstrong
    static constexpr bool enable_3b = false;
    static constexpr bool enable_3b_direct = true;
    static constexpr bool enable_sensor_lj = false;
    static constexpr bool enable_sensor_fene = false;
    static constexpr bool enable_sensor_rdf = false;
    static constexpr bool enable_sensor_disp = false;
    static constexpr bool storeCheckpoint = false;
    static constexpr uint64_t timesteps = 1;            // unit-less
    static constexpr uint64_t write_freq = 10;          // unit-less
    static constexpr int pair_steps = 4;                // unit-less
    static constexpr int triple_steps = 8;              // unit-less

    static constexpr uint64_t sensor_lj_bins = 100;
    static constexpr uint64_t sensor_fene_bins = 100;
    static constexpr bool sensor_temp_opt = false;
    static constexpr double sensor_rdf_max = 5.0;
    static constexpr double sensor_rdf_dr = 0.1;

    static constexpr math::d3 ADR_low {0, 0, 0};
    static constexpr math::d3 ADR_high {0, 0, 0};
    static constexpr math::d3 ADR_h_dim = {0, 0, 0};

    static constexpr uint64_t IBI_bins = 100;
    static constexpr double IBI_alpha = 0.2;
    static constexpr int IBI_steps_equil = 1e+6;
    static constexpr int IBI_steps_measure = 1e+5;
    static constexpr double IBI_conv_threshold = 0.99;
    static constexpr std::string IBI_conv_mode = "integral";
    static constexpr std::string IBI_conv_stop = "worse";
    static constexpr int IBI_conv_window = 10;
    static constexpr math::d3 IBI_exclusion_low {0, 0, 0};
    static constexpr math::d3 IBI_exclusion_high {0, 0, 0};
}

#endif //KOMD_DEFAULTS_H
