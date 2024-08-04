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
            if (var == "domain_low") file >> config->domainLow.x() >> config->domainLow.y() >> config->domainLow.z();
            if (var == "domain_high") file >> config->domainHigh.x() >> config->domainHigh.y() >> config->domainHigh.z();
            if (var == "enableSOA") file >> config->enableSOA;
        }
    }
    file.close();

    return true;
}
