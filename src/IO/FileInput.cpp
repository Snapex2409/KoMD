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
            if (var == "density") file >> config->density;
            if (var == "epsilon") file >> config->epsilon;
            if (var == "sigma") file >> config->sigma;
            if (var == "stiffness") file >> config->stiffness_factor;
            if (var == "limit") file >> config->limit_factor;
            if (var == "domain_low") file >> config->domainLow.x() >> config->domainLow.y() >> config->domainLow.z();
            if (var == "domain_high") file >> config->domainHigh.x() >> config->domainHigh.y() >> config->domainHigh.z();
            if (var == "timesteps") file >> config->timesteps;
            if (var == "write_freq") file >> config->write_freq;
            if (var == "enable_sensor_lj") file >> config->enable_sensor_lj;
            if (var == "enable_sensor_fene") file >> config->enable_sensor_fene;
            if (var == "enable_sensor_rdf") file >> config->enable_sensor_rdf;
            if (var == "sensor_lj_bins") file >> config->sensor_lj_bins;
            if (var == "sensor_fene_bins") file >> config->sensor_fene_bins;
            if (var == "sensor_rdf_max") file >> config->sensor_rdf_max;
            if (var == "sensor_rdf_dr") file >> config->sensor_rdf_dr;
            if (var == "store_checkpoint") file >> config->storeCheckpoint;
            if (var == "checkpoint_file") {
                file >> config->checkpoint_file;
                config->loadCheckpoint = true;
            }
        }
    }
    file.close();

    return true;
}
