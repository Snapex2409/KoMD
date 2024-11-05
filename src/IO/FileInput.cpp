//
// Created by alex on 8/4/24.
//

#include "FileInput.h"
#include "Logging.h"
#include "Registry.h"

#include <fstream>

bool FileInput::readFile(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        Log::io->error() << "Could not open the file: " << filename << std::endl;
        return false;
    }

    auto config = Registry::instance->configuration();
    auto& components = Registry::instance->components();

    const int MAX_LINE_LENGTH = 1024;
    std::string var;
    while (!file.eof() && file.good()) {
        file >> var;
        if (var[0] == '#') { /* ignore comment line*/
            file.ignore(MAX_LINE_LENGTH, '\n');
        } else {
            if (var == "temperature") file >> config->temperature;
            if (var == "delta_t") file >> config->delta_t;
            if (var == "cutoff") file >> config->cutoff;
            if (var == "cell_size") file >> config->cell_size;
            if (var == "density") file >> config->density;
            if (var == "stiffness") file >> config->stiffness_factor;
            if (var == "limit") file >> config->limit_factor;
            if (var == "energy_3b") file >> config->energy_3b;
            if (var == "pair_steps") file >> config->pair_steps;
            if (var == "triple_steps") file >> config->triple_steps;
            if (var == "domain_low") file >> config->domainLow.x() >> config->domainLow.y() >> config->domainLow.z();
            if (var == "domain_high") file >> config->domainHigh.x() >> config->domainHigh.y() >> config->domainHigh.z();
            if (var == "timesteps") file >> config->timesteps;
            if (var == "write_freq") file >> config->write_freq;
            if (var == "enable_3b") file >> config->enable_3b;
            if (var == "enable_sensor_lj") file >> config->enable_sensor_lj;
            if (var == "enable_sensor_fene") file >> config->enable_sensor_fene;
            if (var == "enable_sensor_rdf") file >> config->enable_sensor_rdf;
            if (var == "enable_sensor_disp") file >> config->enable_sensor_disp;
            if (var == "sensor_lj_bins") file >> config->sensor_lj_bins;
            if (var == "sensor_fene_bins") file >> config->sensor_fene_bins;
            if (var == "sensor_rdf_max") file >> config->sensor_rdf_max;
            if (var == "sensor_rdf_dr") file >> config->sensor_rdf_dr;
            if (var == "store_checkpoint") file >> config->storeCheckpoint;
            if (var == "ps_gen_region") {
                math::d3 begin, end; int cid;
                file >> begin.x() >> begin.y() >> begin.z() >> end.x() >> end.y() >> end.z() >> cid;
                if (cid < 0) { Log::io->error() << "negative component id not allowed" << std::endl; return false; }
                config->phasespace_gen_regions.emplace_back(begin, end, cid);
            }
            if (var == "ADR_enable") file >> config->ADR_enable;
            if (var == "ADR_low") file >> config->ADR_low.x() >> config->ADR_low.y() >> config->ADR_low.z();
            if (var == "ADR_high") file >> config->ADR_high.x() >> config->ADR_high.y() >> config->ADR_high.z();
            if (var == "ADR_h_dim") file >> config->ADR_h_dim.x() >> config->ADR_h_dim.y() >> config->ADR_h_dim.z();

            if (var == "IBI_enable") file >> config->IBI_enable;
            if (var == "IBI_bins") file >> config->IBI_bins;
            if (var == "IBI_alpha") file >> config->IBI_alpha;
            if (var == "IBI_steps_equil") {double tmp; file >> tmp; config->IBI_steps_equil; }
            if (var == "IBI_steps_measure") {double tmp; file >> tmp; config->IBI_steps_measure; }
            if (var == "IBI_conv_threshold") file >> config->IBI_conv_threshold;
            if (var == "IBI_conv_mode") file >> config->IBI_conv_mode;
            if (var == "IBI_conv_stop") file >> config->IBI_conv_stop;
            if (var == "IBI_conv_window") file >> config->IBI_conv_window;

            if (var == "checkpoint_file") {
                std::string path; math::d3 offset;
                file >> path >> offset.x() >> offset.y() >> offset.z();
                config->checkpoint_files.emplace_back(path, offset);
            }
            if (var == "enable_one_cell") file >> config->enable_one_cell;

            // Handle component loading
            if (var == "COMP") {
                int cid; double eps, sig, mass; math::d3 r;
                file >> cid >> eps >> sig >> mass >> r.x() >> r.y() >> r.z();
                if (cid < 0) { Log::io->error() << "negative component id not allowed" << std::endl; return false; }
                if (components.size() <= cid) components.resize(cid + 1);
                components[cid].addSite(eps, sig, mass, r);
                if (eps > config->max_epsilon) config->max_epsilon = eps;
                if (sig > config->max_sigma) config->max_sigma = sig;
            }
        }
    }
    file.close();

    // perform some final checks
    if (components.empty()) { Log::io->error() << "no components defined" << std::endl; return false; }
    if (config->max_epsilon == std::numeric_limits<double>::min()) { Log::io->error() << "invalid epsilon specified" << std::endl; return false; }
    if (config->max_sigma == std::numeric_limits<double>::min()) { Log::io->error() << "invalid sigma specified" << std::endl; return false; }
    if (config->phasespace_gen_regions.empty() && config->checkpoint_files.empty()) { Log::io->error() << "No phasespace data specified" << std::endl; return false; }
    for (auto& ps_region : config->phasespace_gen_regions) {
        if (std::get<2>(ps_region) >= components.size()) {
            Log::io->error() << "phasespace gen region uses unknown component" << std::endl;
            return false;
        }
    }
    if (config->cell_size < config->cutoff) {
        Log::io->error() << "cell size must be larger or equal to cutoff radius" << std::endl;
        return false;
    }
    return true;
}
